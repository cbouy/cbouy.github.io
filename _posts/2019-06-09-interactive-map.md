---
image: /assets/blog/2019/06/09/bokeh_plot.png
tags: python bokeh interactive
title: Interactive Choropleth Maps with Bokeh
---
This is an advanced guide to make interactive and animated choropleth maps with Python and the Bokeh library. If you are not familiar with bokeh, please start with their simple tutorials [here](https://bokeh.pydata.org/en/latest/ "Bokeh tutorials").

To recreate the map you will need the following python packages:
* pandas
* geopandas
* bokeh
* matplotlib

The resulting map is a standalone HTML file that you can display on your website or during presentations:

<object data="/assets/blog/2019/06/09/obesity-trends.html" type="text/html" style="overflow:hidden; height: 840px; width: 110%; transform: translateX(-22px);">
    Your browser doesn’t support the object tag.
</object>
To download it: `Right-click > Save link as` [here](/assets/blog/2019/06/09/obesity-trends.html)

You can download my jupyter-notebook with instructions [here](/assets/blog/2019/06/09/obesity-trends.ipynb "Download link").

<object data="/assets/blog/2019/06/09/notebook.html" type="text/html" style="overflow:hidden; height: 800px; width: 100%">
    Your browser doesn’t support the object tag.
</object>

Of course, you can use the same code for any other endpoint. I made a simple reusable function for this in [here](/assets/blog/2019/06/09/interactive_maps.py). Just make sure the dataframe has the same format as the one we used so far, and the column with the numerical value is called "Value" (instead of "Prevalence"). [Here](/assets/blog/2019/06/09/alcohol-consumption.html) is an example using the alcohol consumption per capita.

Cedric
