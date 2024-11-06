"""
Get a table of plasmids queried by "gateway" on AddGene, it takes a table file as input
and appends new plasmids to it.
"""

import httpx
from bs4 import BeautifulSoup
import re
import asyncio


async def main(sort_by, table_file):
    # Read existing plasmids
    existing_plasmids = set()
    with open(table_file, "r") as input_file:
        for line in input_file:
            ls = line.strip().split("\t")
            existing_plasmids.add(ls[0])

    with open(table_file, "a") as output_file:
        for i in range(1, 240):
            print(f"Processing page {i}")
            url = f"https://www.addgene.org/search/catalog/plasmids/?q=gateway&page_number={i}&page_size=50&sort_by={sort_by}"
            async with httpx.AsyncClient() as client:
                resp = await client.get(url)

            soup = BeautifulSoup(resp.content, "html.parser")

            # For each div.search-result-item
            for item in soup.find_all("div", class_="search-result-item"):
                # Extract title from <a href="/227716/">pDest-pol2-Smarcd1-Myc</a>
                a_tag = item.select_one(".search-result-title a")
                if not a_tag:
                    continue
                title = a_tag.text.strip()
                id = re.search(r"/(\d+)/", a_tag["href"])
                if not id:
                    continue
                id = id.group(1)
                if id in existing_plasmids:
                    print(f"Skipping {id} because it already exists")
                    continue
                # Get all the hrefs that start with /browse/article
                article_links = [
                    link["href"].split("/")[-2]
                    for link in item.select("a")
                    if link["href"].startswith("/browse/article")
                ]
                output_file.write(
                    "\t".join([id, title, ",".join(article_links)]) + "\n"
                )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sort_by", type=str, default="newest")
    parser.add_argument(
        "--table-file", type=str, default="data/all_gateway_plasmids.tsv"
    )
    args = parser.parse_args()
    asyncio.run(main(args.sort_by, args.table_file))
