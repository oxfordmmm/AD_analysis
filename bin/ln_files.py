#!/usr/bin/env python3
import os
import sys


def create_symlink(source, link_name):
    try:
        os.symlink(source, link_name)
        print(f"Created symlink: {link_name} -> {source}")
    except FileExistsError:
        print(f"Symlink already exists: {link_name}")
    except OSError as e:
        print(f"Error creating symlink {link_name} -> {source}: {e}")


def main():
    dataFol=sys.argv[1]
    linkFol=sys.argv[2]

    if not os.path.exists(linkFol):
        os.makedirs(linkFol)
    
    for filename in os.listdir(dataFol):
        source_path = os.path.join(dataFol, filename)
        # remove anything before and including '_SQK-RBK114-96_' from filename
        fileRename = filename.split('_SQK-RBK114-96_')[-1]
        #fileRename = filename.replace('7679b82565e0a61142f6119feb7a1f6ba07d4128_SQK-RBK114-96_', '')
        if 'unclassified' in fileRename:
            fileRename = 'barcode00.fastq.gz'
        link_path = os.path.join(linkFol, fileRename)
        create_symlink(source_path, link_path)

if __name__ == "__main__":
    main()