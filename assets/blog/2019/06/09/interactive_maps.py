#!/usr/bin/env python
# coding: utf-8

import textwrap
import pandas as pd
import geopandas as gpd
from bokeh.io import show, output_file, output_notebook, reset_output
from bokeh.plotting import figure
from bokeh.models import (
    GeoJSONDataSource, ColumnDataSource, ColorBar, Slider, Spacer,
    HoverTool, TapTool, Panel, Tabs, Legend, Toggle, LegendItem,
)
from bokeh.palettes import brewer
from bokeh.models.callbacks import CustomJS
from bokeh.models.widgets import Div
from bokeh.layouts import widgetbox, row, column
from matplotlib import pyplot as plt
from matplotlib.colors import rgb2hex

def choropleth_map(
    data, geo,
    bins=[0,2,5,10,15,20,25,30,100], bin_labels="%",
    size=(1100,550), hover_color='#3bdd9d', nodata_color='#d9d9d9',
    output_name="obesity-trends.html", output_title="Obesity trends", notebook=False,
    value_title="Obesity rate (%)", value_axis_label="% Obesity (BMI ≥ 30)",
    map_title="Share of adults who are obese in {}", default_map_year=1975, map_palette='YlOrRd',
    chart_title="Evolution of obesity", default_chart_country="France", country_palette='gist_ncar',
    footer_text=textwrap.dedent("""\
    Data: World Health Organization - Global Health Observatory</br >
    Author: <a href="https://cbouy.github.io">Cédric Bouysset</a>"""),
    ):
    '''
    data: Pandas DataFrame containing the numeric value to display, and the following columns: ['Country', 'Code', 'Year']
    geo: GeoPandas DataFrame containing the info to plot each country. Columns must be name ['Country', 'Code', 'geometry']
    bins: list of bins to label the data. Must contain the min and max possible values.
    bin_labels: labels shown on the graph. Can be a str which will be appended to every bin, or a list of custom labels.
    map_palette: color palette assigned to the bins and displayed on the map. Can be the name of a brewer palette from Bokeh or a list of hex color codes of length len(bins)-1.
    country_palette: color palette used for the chart. Can be the name of a matplotlib colormap or a list of hex color codes.
    '''

    # assign nan values when a year is missing for a country
    df = data[["Year","Code"]].copy()
    mux = pd.MultiIndex.from_product(
        [df["Year"].unique(), df["Code"].unique()],
        names=["Year","Code"]
    )
    df = df.reindex(mux).drop(columns=["Year","Code"]).reset_index()
    data = pd.merge(df, data, how="left", on=["Year","Code"])

    # create stylish labels
    if isinstance(bin_labels, str):
        bin_labels = [f'≤{bins[1]}{bin_labels}'] + [f'{bins[i]}-{bins[i+1]}{bin_labels}' for i in range(1,len(bins)-2)] + [f'>{bins[-2]}{bin_labels}']
    # assign each row to a bin
    data['bin'] = pd.cut(
        data['Value'], bins=bins, right=True, include_lowest=True, precision=0, labels=bin_labels,
    ).astype(str)

    # Merge the geographic data with obesity data
    df = geo.merge(data, on='Code', how='left')
    df = df.drop(columns="Country_y").rename(columns={"Country_x":"Country"})

    # Add a 'No data' bin for countries without data on the endpoint
    df.loc[df['Value'].isna(), 'bin'] = 'No data'
    df.fillna('No data', inplace = True)

    if isinstance(map_palette, str):
        map_palette = brewer[map_palette][len(bins)-1]

    # Assign obesity prevalence to a color
    def val_to_color(value, nan_color=nodata_color):
        if isinstance(value, str): return nan_color
        for i in range(1,len(bins)):
            if value <= bins[i]:
                return map_palette[i-1]
    df['color'] = df['Value'].apply(val_to_color)

    # assign x coordinates
    def bin_to_cbar_x(value):
        if value == 'No data': return -2
        for i,b in enumerate(bin_labels):
            if value == b:
                return 5*(i+1)
    df['cbar_x'] = df['bin'].apply(bin_to_cbar_x)
    # assign width
    df['cbar_w'] = df['Value'].apply(lambda x: 5 if x == 'No data' else 4.7)

    # create color palette for the graph
    countries = sorted(df[df["bin"] != "No data"]["Country"].unique())
    n_country = len(countries)
    if isinstance(country_palette, str):
        cmap = plt.get_cmap(country_palette, n_country)
        country_palette = [rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]

    # define the output
    reset_output()
    if notebook:
        output_notebook(hide_banner=True)
    else:
        output_file(output_name, title=output_title, mode="inline")

    # Input sources
    df.sort_values(by=["Country","Year"], inplace=True)
    # source that will contain all necessary data for the map
    geosource = GeoJSONDataSource(geojson=df.to_json())
    # source that contains the data that is actually shown on the map (for a given year)
    displayed_src = GeoJSONDataSource(geojson=df[df['Year'].isin(['No data', default_map_year])].to_json())
    # source that will be used for the graph (we don't need the countries shapes for this)
    country_source = ColumnDataSource(df[df['Country'] == default_chart_country].drop(columns=["geometry"]))

    # Tools
    # slider to select the year
    temp = pd.to_numeric(df["Year"], errors="coerce")
    min_year = int(temp.min())
    max_year = int(temp.max())
    del temp
    slider = Slider(title='Year',start=min_year, end=max_year, step=1, value=default_map_year)

    # hover tool for the map
    map_hover = HoverTool(tooltips=[
        ('Country','@Country (@Code)'),
        (value_title, '@Value')
    ])

    # hover tool for the graph
    graph_hover = HoverTool(tooltips=[
        ('Country','@Country (@Code)'),
        (value_title, '@Value'),
        ('Year', '@Year')
    ])

    # button for the animation
    anim_button = Toggle(label="▶ Play", button_type="success", width=50, active=False)

    # create map figure
    p = figure(
        title = map_title.format(default_map_year),
        plot_width=size[0], plot_height=size[1],
        toolbar_location="right", tools="tap,pan,wheel_zoom,box_zoom,save,reset", toolbar_sticky=False,
        active_scroll="wheel_zoom",
    )
    p.title.text_font_size = '16pt'
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.axis.visible = False

    # Add hover tool
    p.add_tools(map_hover)

    # Add patches (countries) to the figure
    patches = p.patches(
        'xs','ys', source=displayed_src,
        fill_color='color',
        line_color='black', line_width=0.25, fill_alpha=1,
        hover_fill_color='color',
    )
    # outline when we hover over a country
    patches.hover_glyph.line_color = hover_color
    patches.hover_glyph.line_width = 3
    patches.nonselection_glyph = None

    # create the interactive colorbar
    p_bar = figure(
        title=None, plot_height=80 , plot_width=int(size[0]*0.6),
        tools="tap", toolbar_location=None
    )
    p_bar.xgrid.grid_line_color = None
    p_bar.ygrid.grid_line_color = None
    p_bar.outline_line_color = None
    p_bar.yaxis.visible = False

    # set the title and ticks of the colorbar
    p_bar.xaxis.axis_label = value_axis_label
    p_bar.xaxis.ticker = sorted(df['cbar_x'].unique())
    p_bar.xaxis.major_label_overrides = dict([(i[0],i[1]) for i in df.groupby(['cbar_x','bin']).describe().index])
    p_bar.xaxis.axis_label_text_font_size = "12pt"
    p_bar.xaxis.major_label_text_font_size = "10pt"

    # activate the hover but hide tooltips
    hover_bar = HoverTool(tooltips=None)
    p_bar.add_tools(hover_bar)

    # plot the rectangles for the colorbar
    cbar = p_bar.rect(x='cbar_x', y=0, width='cbar_w', height=1,
        color='color', source=displayed_src,
        hover_line_color=hover_color, hover_fill_color='color')

    # outline when we hover over the colorbar legend
    cbar.hover_glyph.line_width = 4
    cbar.nonselection_glyph = None

    # create the graph figure
    p_country = figure(
        title=chart_title, plot_width=size[0], plot_height=int(size[1]*1.3),
        tools="pan,wheel_zoom,save", active_scroll="wheel_zoom", toolbar_location="right",
    )
    p_country.title.text_font_size = '14pt'
    p_country.xaxis.axis_label = "Year"
    p_country.yaxis.axis_label = value_axis_label
    p_country.axis.major_label_text_font_size = "12pt"
    p_country.axis.axis_label_text_font_size = "14pt"

    # plot data on the figure
    line_plots = {}
    legend_items = {}
    for i, country in enumerate(countries):
        # get subset of data corresponding to a country
        country_source = ColumnDataSource(df[df['Country'] == country].drop(columns=["geometry"]))
        # plot
        line = p_country.line("Year", "Value", legend=False, source=country_source,
                          color=country_palette[i], line_width=2)
        circle = p_country.circle("Year", "Value", legend=False, source=country_source,
                              line_color="darkgrey", fill_color=country_palette[i], size=8)
        # used later in the interactive callbacks
        line_plots[country] = [line, circle]
        legend_items[country] = LegendItem(label=country, renderers=[line, circle])
        # only display France at first
        if country != default_chart_country:
            line.visible = False
            circle.visible = False

    default_legend = [
        (default_chart_country, line_plots[default_chart_country]),
    ]
    legend = Legend(items=default_legend, location="top_center")
    legend.click_policy = "hide"
    p_country.add_layout(legend, 'right')

    # Add hover tool
    p_country.add_tools(graph_hover)

    # JS callbacks

    # Update the map on slider change
    slider_callback = CustomJS(args=dict(slider=slider, source=geosource, displayed_src=displayed_src), code="""
        var year = slider.value;
        var show = [year, 'No data'];
        var data = {};
        columns = Object.keys(source.data);
        columns.forEach(function(key) {
            data[key] = [];
        });
        for (var i = 0; i < source.get_length(); i++){
            if (show.includes(source.data['Year'][i])){
                columns.forEach(function(key) {
                    data[key].push(source.data[key][i])
                });
            }
        }
        displayed_src.data = data;
        displayed_src.change.emit();
    """)
    slider.js_on_change('value', slider_callback)

    # Update figure title from slider change
    callback_title = CustomJS(args=dict(slider=slider, figure=p, title=map_title), code="""
        var year = slider.value;
        figure.title.text = title.replace("{}", year);
    """)
    slider.js_on_change('value', callback_title)


    # Add callback on country click
    plot_callback = CustomJS(args=dict(
        csource=country_source, source=geosource, displayed_src=displayed_src, line_plots=line_plots, legend=legend, legend_items=legend_items), code="""
        // only continue if a country was selected
        var ixs = displayed_src.selected.indices;
        if (ixs.length == 0) { return; }

        // init
        var data = {};
        var items = [];
        countries = [];
        columns = Object.keys(source.data);
        columns.forEach(function(key) {
            data[key] = [];
        });

        // hide all plots
        for (var country in line_plots) {
            var line = line_plots[country][0];
            var circle = line_plots[country][1];
            line.visible = false;
            circle.visible = false;
        }

        // loop over the selected countries
        ixs.forEach(function(ix) {
            // identify corresponding country
            country = displayed_src.data["Country"][ix];
            countries.push(country);
        });
        // sort them in order
        countries.sort()
        // display the corresponding glyphs and legend
        countries.forEach(function(country) {
            line = line_plots[country][0];
            circle = line_plots[country][1];
            line.visible = true;
            circle.visible = true;
            items.push(legend_items[country]);

            for (var i = 0; i < source.get_length(); i++){
                if (source.data['Country'][i] == country) {
                    columns.forEach(function(key) {
                        data[key].push(source.data[key][i])
                    });
                }
            }
        });
        legend.items = items;
        csource.data = data;
        csource.change.emit();
    """)
    displayed_src.selected.js_on_change('indices', plot_callback)

    # add animation
    update_interval = 500 # in ms
    anim_callback = CustomJS(args=dict(slider=slider, update_interval=update_interval, max_year=max_year, min_year=min_year), code="""
        var button = cb_obj;
        if (button.active == true){
            button.label = "◼ Stop";
            button.button_type = "danger";
            mytimer = setInterval(update_year, update_interval);
        } else {
            button.label = "▶ Play";
            button.button_type = "success";
            clearInterval(mytimer);
        }

        function update_year() {
            year = slider.value;
            if (year < max_year) {
                slider.value += 1;
            } else {
                slider.value = min_year;
            }
        }
    """)
    anim_button.callback = anim_callback

    # arrange display with tabs
    tab_map = Panel(title="Map",
        child=column(
            p, # map
            p_bar, # colorbar
            row(widgetbox(anim_button), Spacer(width=10), widgetbox(slider)) # animation button and slider
        ))
    tab_chart = Panel(title="Chart", child=column(p_country))
    tabs = Tabs(tabs=[ tab_map, tab_chart ])

    # save the document and display it !
    footer = Div(text=footer_text)
    layout = column(tabs, footer)
    show(layout)
