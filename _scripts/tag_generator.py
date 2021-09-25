#!/usr/bin/env python

'''
tag_generator.py

Copyright 2017 Long Qian
Contact: lqian8@jhu.edu

Modified by Cédric Bouysset (@cbouy) - 2021

This script creates tags for your Jekyll blog hosted by Github page.
No plugins required.
'''

import os
from pathlib import Path

root_dir = Path(".")
post_dir = root_dir / '_posts'
tag_dir = root_dir / 'tags'

total_tags = []
for filename in post_dir.glob('*md'):
    with open(filename, 'r', encoding='utf-8') as f:
        crawl = False
        for line in f:
            if crawl:
                current_tags = line.strip().split(":")
                if current_tags[0] == 'tags':
                    tags = "".join(current_tags[1:]).split()
                    total_tags.extend(tags)
                    crawl = False
                    break
            if line.strip() == '---':
                if not crawl:
                    crawl = True
                else:
                    crawl = False
                    break
total_tags = set(total_tags)
print(f"{len(total_tags)} unique tags to generate:", total_tags)

if not tag_dir.is_dir():
    tag_dir.mkdir()
else:
    for tag in tag_dir.glob('*.md'):
        os.remove(tag)

for tag in total_tags:
    filename = tag_dir / f"{tag}.md"
    with open(filename, 'w') as f:
        write_str = f'''\
---
layout: tagpage
title: "#{tag}"
tag: {tag}
robots: noindex
---
'''
        f.write(write_str)

n_tags_created = len(list(tag_dir.glob('*.md')))
print(f"Wrote {n_tags_created} tag pages")