import httpx
from bs4 import BeautifulSoup
import re
import asyncio


async def main(sort_by: str = "newest"):
    # Read existing plasmids
    existing_plasmids = set()
    with open("data/all_gateway_plasmids.tsv", "r") as input_file:
        for line in input_file:
            ls = line.strip().split("\t")
            existing_plasmids.add(ls[0])

    with open("data/all_gateway_plasmids.tsv", "a") as output_file:
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
    asyncio.run(main("alpha_desc"))
