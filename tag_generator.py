#!/usr/bin/env python

'''
tag_generator.py

Copyright 2017 Long Qian
Contact: lqian8@jhu.edu

This script creates tags for your Jekyll blog hosted by Github page.
No plugins required.
'''

import glob, os, urllib.parse

post_dir = '_posts/'
tag_dir = 'tags/'

filenames = glob.glob(post_dir + '*md')

total_tags = []
for filename in filenames:
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
print(total_tags)

old_tags = glob.glob(tag_dir + '*.md')
for tag in old_tags:
    os.remove(tag)

if not os.path.exists(tag_dir):
    os.makedirs(tag_dir)

for tag in total_tags:
    tag_filename = tag_dir + tag + '.md'
    with open(tag_filename, 'a') as f:
        write_str = f'''\
---
layout: tagpage
title: "#{tag}"
tag: {tag}
robots: noindex
---
'''
        f.write(write_str)
print("Tags generated, count", total_tags.__len__())
