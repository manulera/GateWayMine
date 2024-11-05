import pandas as pd
import asyncio
import httpx
from bs4 import BeautifulSoup


async def main(input_file, output_file):
    df = pd.read_csv(
        input_file,
        sep="\t",
        names=["addgene_id", "plasmid_name", "article_id"],
        dtype={"addgene_id": str, "article_id": str},
    )
    existing_articles = set()
    with open(output_file, "r") as h:
        for line in h:
            ls = line.strip().split("\t")
            existing_articles.add(ls[0])

    # get unique article ids
    article_ids = list(sorted(df["article_id"].dropna().unique(), key=lambda x: int(x)))

    with open(output_file, "a") as h:
        for article_id in article_ids:
            if article_id in existing_articles:
                print(f"{article_id} already exists")
                continue

            url = f"https://www.addgene.org/browse/article/{article_id}/"
            print(article_id)
            async with httpx.AsyncClient() as client:
                response = await client.get(url)

            soup = BeautifulSoup(response.text, "html.parser")

            # Find first h1 element and get its next p sibling (links are there)
            h1 = soup.find("h1")
            if not h1:
                print(f"No h1 found for {article_id}")
                continue
            p = h1.find_next_sibling("p")
            if not p:
                print(f"No p found for {article_id}")
                continue
            # Get all the hrefs in the p tag
            hrefs = list(sorted(a["href"] for a in p.find_all("a")))
            h.write(f"{article_id}\t{'\t'.join(hrefs)}\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_file", type=str, default="data/all_gateway_plasmids.tsv"
    )
    parser.add_argument(
        "--output_file", type=str, default="data/addgene_article_refs.tsv"
    )
    args = parser.parse_args()
    asyncio.run(main(args.input_file, args.output_file))
