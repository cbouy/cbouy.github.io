---
image: /assets/blog/2022/04/08/banner.png
tags: rdkit python jupyter
title: Fine tuning your Jupyter notebook's displays
---

One of the nice things about Jupyter notebooks is that some objects come with a more visual and practical representation when you put them at the last line of a cell. For example, RDKit displays an image of the molecule, and pandas displays a pretty table view of the dataframe.

In this blog post, I'll show you how to create your own display methods for your classes, and we'll even take it a step further with custom representations for built-in types such as `str` or `list`, and how to enable them for all your notebooks.

<script>
function fit_height(id) {
    var iframe = document.getElementById(id);
    iframe.style.height = 0;
    var height = iframe.contentWindow.document.body.scrollHeight + 20 + 'px';
    iframe.style.height = height;
}
</script>
<iframe id="notebook" src="/assets/blog/2022/04/08/notebook.html" frameborder="0" onload="fit_height('notebook');" width="100%"></iframe>

Here's a gist to my own startup profile:
<div style="max-height: 400px; overflow: auto; margin-bottom: 30px;">
    <script src="https://gist.github.com/cbouy/8406658e74890e08ec8d91ece594e1b8.js"></script>
</div>

For the final *wow* factor, I'll just leave this **abomination** here:

<video autoplay muted controls width="720" src="/assets/blog/2022/04/08/quickmath.mp4">
</video>

If you want to reproduce this (*why would you ever want to do that?*), you'll have to use the `text/plain` output which uses a [slightly different syntax](https://ipython.readthedocs.io/en/stable/api/generated/IPython.lib.pretty.html#module-IPython.lib.pretty) from previous examples, but it has the advantage of working in IPython shells in addition to notebooks.


This concludes the blog post, hope you liked it!

Cheers,

Cédric