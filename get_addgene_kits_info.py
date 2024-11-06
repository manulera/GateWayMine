"""
Gather information about the plasmids contained in AddGene kits used for Gateway cloning
"""

import argparse
import json
import re
import asyncio
from playwright.async_api import async_playwright


async def scrape_addgene_kit(url: str) -> list[str]:
    print(f"Scraping {url}")
    async with async_playwright() as p:
        browser = await p.chromium.launch()
        my_page = await browser.new_page()
        print(f"Navigating to {url + '#kit-contents'}")
        # Change the URL to the site you want for your Playwright web scraping project
        await my_page.goto(url + "#kit-contents")
        await my_page.wait_for_load_state("load")
        print("Page loaded")

        title = (
            (await (await my_page.query_selector("h1#kit-title")).inner_text())
            .strip()
            .split("\n")[0]
        )
        print(f"Title: {title}")

        print("Plasmid table found")

        # Sometimes, there is not #kit-contents, it's split into two, so you have to click
        if not await my_page.query_selector("a[href='#kit-contents']"):
            print("#kit-contents not found, clicking on #kit-contents-1")
            # Click on the link with text "Kit Contents"
            await my_page.click("a[href='#kit-contents-1']")
            await my_page.wait_for_load_state("load")

        # Sometimes not all table is loaded
        if await my_page.query_selector("label:has-text('Show ') >> select"):

            # Error if the option "All" is not found
            if not await my_page.query_selector(
                "label:has-text('Show ') >> select >> option[value='-1']"
            ):
                raise Exception(
                    "Problem scraping AddGene website: option 'All' not found in the table of plasmids"
                )

            # Click on a select element inside a label
            # the label contains text "Show "
            # include label in the query
            # Select option "All"
            await my_page.select_option("label:has-text('Show ') >> select", "All")
            await my_page.wait_for_load_state("load")

        # From the table (div.panel-body.inventory-panel, get all links (td span a))

        # Use CSS selector to get all the link elements within the specific div
        links = await my_page.query_selector_all(
            "div.panel-body.inventory-panel table.datatable td span a"
        )
        if len(links) == 0:
            links = await my_page.query_selector_all(".tab-pane.active table td a")

        # Extract and print the href attribute from each link
        hrefs = [await link.get_attribute("href") for link in links]
        names = list()
        addgene_ids = list()
        for i, href in enumerate(hrefs):
            match = re.search(r"/(\d+)/", href)
            if match is not None:
                addgene_ids.append(match.group(1))
                name = await links[i].inner_text()
                # Sometimes the link is on the addgene_id, so we need to get the name from the next td
                if re.match(r"^\d+$", name):
                    # Get parent td of the link, and find next sibling td
                    parent_td = await links[i].evaluate_handle(
                        'element => element.closest("td")'
                    )

                    # Find the next sibling td (or handle as required)
                    next_td = await parent_td.evaluate_handle(
                        "td => td.nextElementSibling"
                    )
                    name = await next_td.inner_text()

                names.append(name)
            else:
                print(f"skipped link {href}")
        print()
        print()
        # Close the browser
        await browser.close()

    return title, addgene_ids, names


async def main(input_file, output_file):
    urls = list()
    with open(input_file, "r") as handle:
        for line in handle:
            urls.append(line.strip())
    out_dict = dict()
    for url in urls:
        title, addgene_ids, names = await scrape_addgene_kit(url)
        out_dict[title] = dict()
        out_dict[title]["name"] = title
        out_dict[title]["url"] = url
        out_dict[title]["addgene_ids"] = addgene_ids
        out_dict[title]["plasmid_names"] = names
    with open(output_file, "w") as handle:
        json.dump(out_dict, handle, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input_file",
        required=False,
        help="File with the URLs of the AddGene kits",
        default="data/addgene_kits.txt",
    )
    parser.add_argument(
        "--output_file",
        required=False,
        help="File to save the plasmid files",
        default="data/addgene_kit_plasmids.json",
    )
    args = parser.parse_args()
    asyncio.run(main(args.input_file, args.output_file))
