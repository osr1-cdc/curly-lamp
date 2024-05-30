library(networkD3)
library(ape)
library(data.tree)
library(tidytree)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(webshots) # Also run: webshot::install_phantomjs()
library(htmlwidgets)


# CSV file with lineage colors
colors_csv <- "Lineage_Colors_on_CDT.csv"

# manual lineage relationship designation
clean_string_edit <- c("Origin Strain/B.1.617.2",
                       "Origin Strain/B.1.1.529",
                       "Origin Strain/B.1.1.529/BA.1/BA.1.1",
                       "Origin Strain/B.1.1.529/BA.2",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.12.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.75",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.75/BA.2.75.2",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.75/CH.1.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.75/BN.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.11.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.11.1/KP.1.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.11.1/KP.2",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.11.1/KP.3",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.13",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.13/JN.1.13.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.13/JN.1.13.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.16",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.16/JN.1.16.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.18",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/KW.1.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.32",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/KQ.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/KV.2",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.7",
                       "Origin Strain/B.1.1.529/BA.2/BA.2.86/JN.1/JN.1.8.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.4",
                       "Origin Strain/B.1.1.529/BA.2/BA.4/BA.4.6",
                       "Origin Strain/B.1.1.529/BA.2/BA.5",
                       "Origin Strain/B.1.1.529/BA.2/BA.5/BF.11",
                       "Origin Strain/B.1.1.529/BA.2/BA.5/BF.7",
                       "Origin Strain/B.1.1.529/BA.2/BA.5/BA.5.2.6",
                       "Origin Strain/B.1.1.529/BA.2/BA.5/BQ.1",
                       "Origin Strain/B.1.1.529/BA.2/BA.5/BQ.1/BQ.1.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.16",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.16/XBB.1.16.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.16/XBB.1.16.11",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.16/HF.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.16/XBB.1.16.15",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.16/XBB.1.16.17",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.16/XBB.1.16.6",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.16/JF.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/FE.1.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.42.2",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/XBB.1.5.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/XBB.1.5.10",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/JD.1.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/FD.1.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/FD.2",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/EU.1.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/XBB.1.5.59",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/XBB.1.5.68",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/XBB.1.5.70",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/XBB.1.5.70/GK.1.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/XBB.1.5.70/GK.2",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.5/XBB.1.5.70/XBB.1.5.72",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.9.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.9.1/FL.1.5.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.9.2",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.9.2/EG.5",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.9.2/EG.5/HK.3",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.9.2/EG.5/JG.3",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.9.2/EG.5/HV.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.9.2/EG.5/EG.5.1.8",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.1.9.2/EG.6.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.2.3",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.2.3/GE.1",
                       "Origin Strain/B.1.1.529/BA.2/XBB/XBB.2.3/XBB.2.3.8",
                       "Origin Strain/B.1.1.529/BA.2/XDP")

clean_string_df_filter <- tibble(clean_string = clean_string_edit) %>%
  mutate(node_name = str_match(clean_string, ".+/([^/]+$)")[,2])

clean_string_df_filter$pathString <- clean_string_df_filter$clean_string

pango_node_filter <- as.Node(clean_string_df_filter)

pagnotreeList <- ToListExplicit(pango_node_filter, unname = TRUE)


# Default coloring list for lineages
force_internal_colors <- FALSE
internal_color_list <- list(`Origin Strain` = "#797979",
                   B.1.1.529	= "#E26028",
                   BA.1	= "#E26028",
                   BA.1.1	= "#FF824C",
                   BA.2	= "#9CCC65",
                   BA.2.12.1	= "#7CB342",
                   BA.2.75	= "#D4E157",
                   BA.2.75.2	= "#C0CA33",
                   CH.1.1	= "#827717",
                   BN.1	= "#9e9d24",
                   BA.4	= "#FFD54F",
                   BA.4.6	= "#FFB300",
                   BA.5	= "#80CBC4",
                   BF.7	= "#81D4FA",
                   BF.11	= "#29B6F6",
                   BA.5.2.6	= "#009688",
                   BQ.1	= "#006064",
                   BQ.1.1	= "#00838F",
                   XBB	= "#9FA8DA",
                   XBB.1.5	= "#3F51B5 ",
                   XBB.1.5.1	= "#1a237e",
                   XBB.1.5.10	= "#f97a57",
                   EU.1.1	= "#f25353",
                   XBB.1.5.68	= "#86bcb6",
                   XBB.1.9.1	= "#304ffe",
                   XBB.1.9.2	= "#536DFE",
                   FD.2	= "#8C9EFF",
                   XBB.1.16	= "#4527A0",
                   XBB.1.16.1	= "#f89b9b",
                   FE.1.1	= "#499894",
                   XBB.2.3	= "#f7bcdf",
                   B.1.617.2	= "#B39DDB",
                   Other	= "#797979",
                   EG.5	= "#f28e2b",
                   XBB.1.5.59	= "#ffbe7d",
                   XBB.1.16.6	= "#59a14f",
                   XBB.1.5.72	= "#254220",
                   FL.1.5.1	= "#8cd17d",
                   GE.1	= "#66534c",
                   EG.6.1	= "#b9aba8",
                   XBB.1.16.11	= "#d37295",
                   FD.1.1	= "#e7cb40",
                   XBB.1.5.70	= "#46783c",
                   XBB.2.3.8	= "#b9b2e8",
                   HV.1	= "#797765",
                   XBB.1.42.2	= "#75d1a3",
                   XBB.1.16.15	= "#e3b5fa",
                   GK.2	= "#9d8577",
                   HF.1	= "#d399ad",
                   GK.1.1	= "#90fae5",
                   HK.3	= "#859dbc",
                   JD.1.1	= "#80f192",
                   JF.1	= "#fff1a9",
                   JG.3	= "#66a172",
                   BA.2.86	= "#d771f1",
                   JN.1	= "#660099",
                   XBB.1.16.17	= "#8f8cb0",
                   EG.5.1.8	= "#e14d6f",
                   JN.1.13	= "#e83c78",
                   JN.1.18	= "#4df230",
                   JN.1.16	= "#bc3a88",
                   JN.1.7	= "#cf15e8",
                   JN.1.8.1	= "#ffb027",
                   KQ.1	= "#e86cde",
                   KP.2	= "#5875a1",
                   JN.1.13.1	= "#f2591c",
                   KP.1.1	= "#cbff86",
                   JN.1.11.1	= "#5c9cf1",
                   XDP	= "#b60052",
                   JN.1.16.1	= "#366298",
                   KP.3	= "#b06589",
                   JN.1.32	= "#2e7919",
                   KW.1.1	= "#34a0ff",
                   KS.1	= "#b771bc",
                   KV.2	= "#e14365"
                   )




# Check if the CSV file exists
if (file.exists(colors_csv) && force_internal_colors != TRUE) {
  # Step 1: Read the CSV file into a data frame
  print(paste0(colors_csv, " found! Using it."))
  color_data <- read.csv(colors_csv, stringsAsFactors = FALSE)
  
  # Step 2: Convert the data frame to a named list
  color_list <- setNames(as.list(color_data$Hex_Value), color_data$Variant)
} else {
  # Use the default list
  color_list <- internal_color_list
  print("Using internal color set")
}



prop_list <-list(HV.1 = 0.289585804,
                 XBB.1.16.11 = 0.025101432,
                 HK.3 = 0.078466785,
                 XBB.1.5 = 0.003074311,
                 XBB = 0.007794122,
                 HF.1 = 0.013844966,
                 XBB.1.9.1 = 0.003044804,
                 JD.1.1 = 0.045540183,
                 GK.1.1 = 0.016412872,
                 GE.1 = 0.008551353,
                 XBB.1.16.15 = 0.013531104,
                 EG.6.1 = 0.004821252,
                 GK.2 = 0.006405419,
                 XBB.1.9.2 = 0.00161383,
                 CH.1.1 = 0.001606371,
                 BA.2 = 0.010081995,
                 XBB.1.5.68 = 0.001975924,
                 XBB.1.5.72 = 0.001607387,
                 XBB.1.42.2 = 0.001987982,
                 XBB.2.3.8 = 0.001088461,
                 XBB.1.5.10 = 0.000659756,
                 XBB.1.5.59 = 0.000328166,
                 FD.1.1 = 0.000274138,
                 EU.1.1 = 4.70E-05,
                 XBB.1.5.1 = 0.000116715,
                 FE.1.1 = 6.00E-05,
                 BA.4.6 = 0.000337744,
                 BA.4 = 0.000335349,
                 BQ.1.1 = 7.16E-07,
                 BA.2.75 = 5.52E-05,
                 BA.1.1 = 2.20E-05,
                 EG.5 = 0.216580613,
                 XBB.2.3 = 0.018538159,
                 XBB.1.16.6 = 0.090533013,
                 FL.1.5.1 = 0.09273139,
                 XBB.1.16.1 = 0.008088904,
                 XBB.1.5.70 = 0.011916313,
                 XBB.1.16 = 0.012365146,
                 Other = 0.010873377
  
)
#prop_list <-list()
# coloring frame and functions
color.df <- tibble(name = names(color_list), color = as.character(color_list))
w <- paste('{', paste(color.df %>% 
                        mutate(name = paste0('"', name), color = paste0(color, '"')) %>%
                        unite('x', c(name, color), sep = '" : "' ) %>%
                        .$x, collapse= ', '), '}', collapse = '')

node.col.func <- JS(paste0('function(d, i) { return ', w, '[d.data.name]; }'))
rm(net_obj)
net_obj <- diagonalNetwork(pagnotreeList, 
                fontFamily = "arial",
                fontSize = 16,
                nodeColour = node.col.func,
                textColour = node.col.func,
                nodeStroke = node.col.func,
                linkColour = "#999",
                opacity = 1.0,
                height= 700, width=900,
                margin=list(left=1,top=0, bottom=0, right=1))

# write some JavaScript to set the circles' radii to the data's value
node_size_js <-
  '
    function(el) { 
      d3.select(el)
        .selectAll(".node")
        .selectAll("circle")
        .attr("r", d => 6);
    }
  '

# use onRender to run the JS when it loads
htmlwidgets::onRender(net_obj, node_size_js)


######

# Create the HTML file for the htmlwidget
html_file <- "pango_tree_diagram.html"
saveWidget(net_obj, file = html_file)

# Capture a screenshot of the htmlwidget and save it as a PNG image
png_file <- "pango_tree_diagram.png"
webshot(url = paste0("file://", normalizePath(html_file)), 
        vwidth = 1000,
        vheight = 700, zoom=2,
        file = png_file)

# Remove the temporary HTML file
file.remove(html_file)

