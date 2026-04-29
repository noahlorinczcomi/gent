#!/usr/bin/env python3
"""Download and extract gent LD reference matrices from Zenodo."""

import argparse
import os
import tarfile
import urllib.request
import zipfile

# EUR/EAS/SAS/AMR are true zip archives; AFR parts are gzip-compressed tarballs
# despite the .zip extension, so they require tar extraction.
URLS = {
    "EUR": [("zip", "https://zenodo.org/records/16762678/files/EUR.zip")],
    "EAS": [("zip", "https://zenodo.org/records/16782543/files/EAS.zip")],
    "SAS": [("zip", "https://zenodo.org/records/16755767/files/SAS.zip")],
    "AMR": [("zip", "https://zenodo.org/records/16782541/files/AMR.zip")],
    "AFR": [
        ("tar", "https://zenodo.org/records/16780125/files/AFR_chr1-10.zip"),
        ("tar", "https://zenodo.org/records/16782180/files/AFR_chr11-22.zip"),
    ],
}


def _progress(block_num, block_size, total_size):
    downloaded = block_num * block_size
    if total_size > 0:
        pct = min(100.0, downloaded * 100 / total_size)
        print(f"\r  {pct:5.1f}%  {downloaded/1e6:.1f} / {total_size/1e6:.1f} MB", end="", flush=True)
    else:
        print(f"\r  {downloaded/1e6:.1f} MB", end="", flush=True)


def fetch(fmt, url, dest_dir):
    filename = os.path.basename(url)
    dest_path = os.path.join(dest_dir, filename)

    print(f"Downloading {filename} ...")
    urllib.request.urlretrieve(url, dest_path, reporthook=_progress)
    print()

    print(f"Extracting {filename} ...")
    if fmt == "tar":
        with tarfile.open(dest_path, "r:gz") as tf:
            tf.extractall(dest_dir)
    else:
        with zipfile.ZipFile(dest_path, "r") as zf:
            zf.extractall(dest_dir)
    os.remove(dest_path)
    print(f"Done: {filename}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Download gent LD reference matrices from Zenodo.",
        usage="pixi run fetch-ld <population> [output_dir]",
    )
    parser.add_argument(
        "population",
        metavar="population",
        help="Population to download: EUR, EAS, SAS, AMR, AFR, or ALL (case-insensitive)",
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        default=".",
        help="Directory in which to place the extracted LD matrices (default: current directory)",
    )
    args = parser.parse_args()

    pop = args.population.upper()
    valid = set(URLS) | {"ALL"}
    if pop not in valid:
        parser.error(f"Invalid population '{args.population}'. Choose from: {', '.join(sorted(valid))}.")

    dest = args.output_dir
    os.makedirs(dest, exist_ok=True)

    targets = list(URLS) if pop == "ALL" else [pop]
    for target in targets:
        for fmt, url in URLS[target]:
            fetch(fmt, url, dest)

    print("All downloads complete.")


if __name__ == "__main__":
    main()
