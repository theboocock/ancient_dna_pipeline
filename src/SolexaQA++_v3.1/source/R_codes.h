#ifndef R_CODE
#define R_CODE
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>



int R_heatmap(char* directory, std::string filepath, std::string filename, std::string matrix_filename_str){

    const char* R_file;
    std::string R_filename_str;
    if(directory){
        R_filename_str=directory+filename+".R";
        R_file=R_filename_str.c_str();
    }else{
        R_filename_str=filepath+".R";
        R_file=R_filename_str.c_str();
    }

    std::ofstream R (R_file, std::ofstream::out);
    if(not R.is_open()){
        std::cout<<"Error: Failure opening "<<R_file<<" for writing."<<std::endl;
        return 1;
    }

    for (int i = 0; i < matrix_filename_str.length(); ++i)
        if (matrix_filename_str[i] == '\\')
        {
        matrix_filename_str.insert(i, 1, '\\');
        ++i; // Skip inserted char
        }

    R<<"library(package=grid)\n"
        "lo = function(rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, treeheight_col, treeheight_row, legend, annotation, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, ...){\n"
        "# Get height of colnames and length of rownames\n"
        "if(!is.null(coln[1])){\n"
        "longest_coln = which.max(nchar(coln))\n"
        "gp = list(fontsize = fontsize_col, ...)\n"
        "coln_height = unit(1.1, \"grobheight\", textGrob(coln[longest_coln], rot = 90, gp = do.call(gpar, gp)))\n"
        "}\n"
        "else{\n"
        "coln_height = unit(5, \"bigpts\")\n"
        "}\n\n"

        "if(!is.null(rown[1])){\n"
        "longest_rown = which.max(nchar(rown))\n"
        "gp = list(fontsize = fontsize_row, ...)\n"
        "rown_width = unit(1.2, \"grobwidth\", textGrob(rown[longest_rown], gp = do.call(gpar, gp)))\n"
        "}\n"
        "else{\n"
        "rown_width = unit(5, \"bigpts\")\n"
        "}\n"

        "gp = list(fontsize = fontsize, ...)\n"
        "# Legend position\n"
        "if(!is.na(legend[1])){\n"
        "longest_break = which.max(nchar(names(legend)))\n"
        "longest_break = unit(1.1, \"grobwidth\", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))\n"
        "title_length = unit(1.1, \"grobwidth\", textGrob(\"Scale\", gp = gpar(fontface = \"bold\", ...)))\n"
        "legend_width = unit(12, \"bigpts\") + longest_break * 1.2\n"
        "legend_width = max(title_length, legend_width)\n"
        "}\n"
        "else{\n"
        "legend_width = unit(0, \"bigpts\")\n"
        "}\n"
        "\n"
        "# Set main title height\n"
        "if(is.na(main)){\n"
        "main_height = unit(0, \"npc\")\n"
        "}\n"
        "else{\n"
        "main_height = unit(5.5, \"grobheight\", textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))\n"
        "}\n"
        "\n"
        "# Column annotations\n"
        "if(!is.na(annotation[[1]][1])){\n"
        "# Column annotation height\n"
        "annot_height = unit(ncol(annotation) * (8 + 2) + 2, \"bigpts\")\n"
        "# Width of the correponding legend\n"
        "longest_ann = which.max(nchar(as.matrix(annotation)))\n"
        "annot_legend_width = unit(1.2, \"grobwidth\", textGrob(as.matrix(annotation)[longest_ann], gp = gpar(...))) + unit(12, \"bigpts\")\n"
        "if(!annotation_legend){\n"
        "annot_legend_width = unit(0, \"npc\")\n"
        "}\n"
        "}\n"
        "else{\n"
        "annot_height = unit(0, \"bigpts\")\n"
        "annot_legend_width = unit(0, \"bigpts\")\n"
        "}\n"
        "\n"
        "# Tree height\n"
        "treeheight_col = unit(treeheight_col, \"bigpts\") + unit(5, \"bigpts\")\n"
        "treeheight_row = unit(treeheight_row, \"bigpts\") + unit(5, \"bigpts\")\n"
        "\n"
        "# Set cell sizes\n"
        "if(is.na(cellwidth)){\n"
        "matwidth = unit(1, \"npc\") - rown_width - legend_width - treeheight_row - annot_legend_width\n"
        "}\n"
        "else{\n"
        "matwidth = unit(cellwidth * ncol, \"bigpts\")\n"
        "}\n"
        "\n"
        "if(is.na(cellheight)){\n"
        "matheight = unit(1, \"npc\") - main_height - coln_height - treeheight_col - annot_height\n"
        "}\n"
        "else{\n"
        "matheight = unit(cellheight * nrow, \"bigpts\")\n"
        "}\n"
        "\n"
        "\n"
        "# Produce layout()\n"
        "pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 6, widths = unit.c(treeheight_row, matwidth, rown_width, legend_width, annot_legend_width), heights = unit.c(main_height, treeheight_col, annot_height, matheight, coln_height)), gp = do.call(gpar, gp)))\n"
        "\n"
        "# Get cell dimensions\n"
        "pushViewport(vplayout(4, 2))\n"
        "cellwidth = convertWidth(unit(0:1, \"npc\"), \"bigpts\", valueOnly = T)[2] / ncol\n"
        "cellheight = convertHeight(unit(0:1, \"npc\"), \"bigpts\", valueOnly = T)[2] / nrow\n"
        "upViewport()\n"
        "\n"
        "# Return minimal cell dimension in bigpts to decide if borders are drawn\n"
        "mindim = min(cellwidth, cellheight)\n"
        "return(mindim)\n"
        "}\n"
        "\n"
        "draw_dendrogram = function(hc, horizontal = T){\n"
        "h = hc$height / max(hc$height) / 1.05\n"
        "m = hc$merge\n"
        "o = hc$order\n"
        "n = length(o)\n"
        "\n"
        "m[m > 0] = n + m[m > 0]\n"
        "m[m < 0] = abs(m[m < 0])\n"
        "\n"
        "dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c(\"x\", \"y\")))\n"
        "dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)\n"
        "\n"
        "for(i in 1:nrow(m)){\n"
        "dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2\n"
        "dist[n + i, 2] = h[i]\n"
        "}\n"
        "\n"
        "draw_connection = function(x1, x2, y1, y2, y){\n"
        "grid.lines(x = c(x1, x1), y = c(y1, y))\n"
        "grid.lines(x = c(x2, x2), y = c(y2, y))\n"
        "grid.lines(x = c(x1, x2), y = c(y, y))\n"
        "}\n"
        "\n"
        "if(horizontal){\n"
        "for(i in 1:nrow(m)){\n"
        "draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])\n"
        "}\n"
        "}\n"
        "\n"
        "else{\n"
        "gr = rectGrob()\n"
        "pushViewport(viewport(height = unit(1, \"grobwidth\", gr), width = unit(1, \"grobheight\", gr), angle = 90))\n"
        "dist[, 1] = 1 - dist[, 1]\n"
        "for(i in 1:nrow(m)){\n"
        "draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])\n"
        "}\n"
        "upViewport()\n"
        "}\n"
        "}\n"
        "\n"
        "draw_matrix = function(matrix, border_color, fmat, fontsize_number){\n"
        "n = nrow(matrix)\n"
        "m = ncol(matrix)\n"
        "x = (1:m)/m - 1/2/m\n"
        "y = 1 - ((1:n)/n - 1/2/n)\n"
        "for(i in 1:m){\n"
        "grid.rect(x = x[i], y = y[1:n], width = 1/m, height = 1/n, gp = gpar(fill = matrix[,i], col = border_color))\n"
        "if(attr(fmat, \"draw\")){\n"
        "grid.text(x = x[i], y = y[1:n], label = fmat[, i], gp = gpar(col = \"grey30\", fontsize = fontsize_number))\n"
        "}\n"
        "}\n"
        "}\n"
        "\n"
        "draw_colnames = function(coln, ...){\n"
        "m = length(coln)\n"
        "x = (1:m)/m - 1/2/m\n"
        "grid.text(coln, x = x, y = unit(0.96, \"npc\"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))\n"
        "}\n"
        "\n"
        "draw_rownames = function(rown, ...){\n"
        "n = length(rown)\n"
        "y = 1 - ((1:n)/n - 1/2/n)\n"
        "grid.text(rown, x = unit(0.04, \"npc\"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))\n"
        "}\n"
        "draw_titles = function(rown, ...){\n"
        "n = length(rown)\n"
        "y = 3 - ((1:n)/n - 1/2/n)\n"
        "grid.text(rown, x = unit(0.04, \"npc\"), y = y, vjust = 0.5, hjust = 0, gp = gpar(fontface='bold', ...))\n"
        "}\n"
        "\n"
        "draw_legend = function(color, breaks, legend, ...){\n"
        "height = min(unit(1, \"npc\"), unit(150, \"bigpts\"))\n"
        "pushViewport(viewport(x = 0, y = unit(1, \"npc\"), just = c(0, 1), height = height))\n"
        "legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))\n"
        "breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))\n"
        "h = breaks[-1] - breaks[-length(breaks)]\n"
        "grid.rect(x = 0, y = breaks[-length(breaks)], width = unit(10, \"bigpts\"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = \"#FFFFFF00\"))\n"
        "grid.text(names(legend), x = unit(12, \"bigpts\"), y = legend_pos, hjust = 0, gp = gpar(...))\n"
        "upViewport()\n"
        "}\n"
        "\n"
        "convert_annotations = function(annotation, annotation_colors){\n"
        "new = annotation\n"
        "for(i in 1:ncol(annotation)){\n"
        "a = annotation[, i]\n"
        "b = annotation_colors[[colnames(annotation)[i]]]\n"
        "if(is.character(a) | is.factor(a)){\n"
        "a = as.character(a)\n"
        "if(length(setdiff(a, names(b))) > 0){\n"
        "stop(sprintf(\"Factor levels on variable %s do not match with annotation_colors\", colnames(annotation)[i]))\n"
        "}\n"
        "new[, i] = b[a]\n"
        "}\n"
        "else{\n"
        "a = cut(a, breaks = 100)\n"
        "new[, i] = colorRampPalette(b)(100)[a]\n"
        "}\n"
        "}\n"
        "return(as.matrix(new))\n"
        "}\n"
        "\n"
        "draw_annotations = function(converted_annotations, border_color){\n"
        "n = ncol(converted_annotations)\n"
        "m = nrow(converted_annotations)\n"
        "x = (1:m)/m - 1/2/m\n"
        "y = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))\n"
        "for(i in 1:m){\n"
        "grid.rect(x = x[i], unit(y[1:n], \"bigpts\"), width = 1/m, height = unit(8, \"bigpts\"), gp = gpar(fill = converted_annotations[i, ], col = border_color))\n"
        "}\n"
        "}\n"
        "\n"
        "draw_annotation_legend = function(annotation, annotation_colors, border_color, ...){\n"
        "y = unit(1, \"npc\")\n"
        "text_height = unit(1, \"grobheight\", textGrob(\"FGH\", gp = gpar(...)))\n"
        "for(i in names(annotation_colors)){\n"
        "grid.text(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = \"bold\", ...))\n"
        "y = y - 1.5 * text_height\n"
        "if(is.character(annotation[, i]) | is.factor(annotation[, i])){\n"
        "for(j in 1:length(annotation_colors[[i]])){\n"
        "grid.rect(x = unit(0, \"npc\"), y = y, hjust = 0, vjust = 1, height = text_height, width = text_height, gp = gpar(col = border_color, fill = annotation_colors[[i]][j]))\n"
        "grid.text(names(annotation_colors[[i]])[j], x = text_height * 1.3, y = y, hjust = 0, vjust = 1, gp = gpar(...))\n"
        "y = y - 1.5 * text_height\n"
        "}\n"
        "}\n"
        "else{\n"
        "yy = y - 4 * text_height + seq(0, 1, 0.02) * 4 * text_height\n"
        "h = 4 * text_height * 0.02\n"
        "grid.rect(x = unit(0, \"npc\"), y = yy, hjust = 0, vjust = 1, height = h, width = text_height, gp = gpar(col = \"#FFFFFF00\", fill = colorRampPalette(annotation_colors[[i]])(50)))\n"
        "txt = rev(range(grid.pretty(range(annotation[, i], na.rm = TRUE))))\n"
        "yy = y - c(0, 3) * text_height\n"
        "grid.text(txt, x = text_height * 1.3, y = yy, hjust = 0, vjust = 1, gp = gpar(...))\n"
        "y = y - 4.5 * text_height\n"
        "}\n"
        "y = y - 1.5 * text_height\n"
        "}\n"
        "}\n"
        "draw_main = function(text, ...){\n"
        "grid.text(text, gp = gpar(fontface = \"bold\", ...))\n"
        "}\n"
        "\n"
        "vplayout = function(x, y){\n"
        "return(viewport(layout.pos.row = x, layout.pos.col = y))\n"
        "}\n"
        "\n"
        "heatmap_motor = function(matrix, border_color, cellwidth, cellheight, tree_col, tree_row, treeheight_col, treeheight_row, filename, width, height, breaks, color, legend, annotation, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, fmat, fontsize_number, ...){\n"
        "grid.newpage()\n"
        "\n"
        "# Set layout\n"
        "mindim = lo(coln = colnames(matrix), rown = rownames(matrix), nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col,  ...)\n"
        "\n"
        "if(!is.na(filename)){\n"
        "pushViewport(vplayout(1:5, 1:5))\n"
        "\n"
        "if(is.na(height)){\n"
        "height = convertHeight(unit(0:1, \"npc\"), \"inches\", valueOnly = T)[2]\n"
        "}\n"
        "if(is.na(width)){\n"
        "width = convertWidth(unit(0:1, \"npc\"), \"inches\", valueOnly = T)[2]\n"
        "}\n"
        "\n"
        "# Get file type\n"
        "r = regexpr(\"\\\\.[a-zA-Z]*$\", filename)\n"
        "if(r == -1) stop(\"Improper filename\")\n"
        "ending = substr(filename, r + 1, r + attr(r, \"match.length\"))\n"
        "\n"
        "f = switch(ending,\n"
        "pdf = function(x, ...) pdf(x, ...),\n"
        "png = function(x, ...) png(x, units = \"in\", res = 300, ...),\n"
        "jpeg = function(x, ...) jpeg(x, units = \"in\", res = 300, ...),\n"
        "jpg = function(x, ...) jpeg(x, units = \"in\", res = 300, ...),\n"
        "tiff = function(x, ...) tiff(x, units = \"in\", res = 300, compression = \"lzw\", ...),\n"
        "bmp = function(x, ...) bmp(x, units = \"in\", res = 300, ...),\n"
        "stop(\"File type should be: pdf, png, bmp, jpg, tiff\")\n"
        ")\n"
        "\n"
        "# print(sprintf(\"height:%f width:%f\", height, width))\n"
        "f(filename, height = height+2, width = width+1)\n"
        "heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, ...)\n"
        "garbage<-dev.off()\n"
        "upViewport()\n"
        "return()\n"
        "}\n"
        "\n"
        "# Omit border color if cell size is too small\n"
        "if(mindim < 3) border_color = NA\n"
        "\n"
        "# Draw title\n"
        "if(!is.na(main)){\n"
        "pushViewport(vplayout(1, 2))\n"
        "draw_main(main, fontsize = 1.6 * fontsize, ...)\n"
        "upViewport()\n"
        "}\n"

        "\n"
        "# Draw tree for the columns\n"
        "if(!is.na(tree_col[[1]][1]) & treeheight_col != 0){\n"
        "pushViewport(vplayout(2, 2))\n"
        "draw_dendrogram(tree_col, horizontal = T)\n"
        "upViewport()\n"
        "}\n"
        "\n"
        "# Draw tree for the rows\n"
        "if(!is.na(tree_row[[1]][1]) & treeheight_row != 0){\n"
        "pushViewport(vplayout(4, 1))\n"
        "draw_dendrogram(tree_row, horizontal = F)\n"
        "upViewport()\n"
        "}\n"
        "\n"
        "# Draw matrix\n"
        "pushViewport(vplayout(4, 2))\n"
        "draw_matrix(matrix, border_color, fmat, fontsize_number)\n"
        "upViewport()\n"
        "\n"
        "# Draw colnames\n"
        "if(length(colnames(matrix)) != 0){\n"
        "pushViewport(vplayout(5, 2))\n"
        "pars = list(colnames(matrix), fontsize = fontsize_col, ...)\n"
        "do.call(draw_colnames, pars)\n"
        "upViewport()\n"
        "}\n"
        "\n"
        "# Draw rownames\n"
        "if(length(rownames(matrix)) != 0){\n"
        "pushViewport(vplayout(4, 3))\n"
        "pars = list(rownames(matrix), fontsize = fontsize_row, ...)\n"
        "do.call(draw_rownames, pars)\n"
        "upViewport()\n"
        "}\n\n"
        "# Draw titles\n"
        "pushViewport(vplayout(2, 3))\n"
        "pars = list(c(\"Tile  (Density)\"), fontsize = fontsize_row, ...)\n"
        "do.call(draw_titles, pars)\n"
        "upViewport()\n"
        "\n"
        "# Draw annotation tracks\n"
        "if(!is.na(annotation[[1]][1])){\n"
        "pushViewport(vplayout(3, 2))\n"
        "converted_annotation = convert_annotations(annotation, annotation_colors)\n"
        "draw_annotations(converted_annotation, border_color)\n"
        "upViewport()\n"
        "}\n"
        "\n"
        "# Draw annotation legend\n"
        "if(!is.na(annotation[[1]][1]) & annotation_legend){\n"
        "if(length(rownames(matrix)) != 0){\n"
        "pushViewport(vplayout(4:5, 5))\n"
        "}\n"
        "else{\n"
        "pushViewport(vplayout(3:5, 5))\n"
        "}\n"
        "draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)\n"
        "upViewport()\n"
        "}\n"
        "\n"
        "# Draw legend\n"
        "if(!is.na(legend[1])){\n"
        "length(colnames(matrix))\n"
        "if(length(rownames(matrix)) != 0 && length(rownames(matrix))>3){\n"
        "pushViewport(vplayout(4:5, 5))\n"
        "}else if(length(rownames(matrix))!= 0 && length(rownames(matrix))<3){\n"
        "pushViewport(vplayout(1:5, 5))\n"
        "}else{\n"
        "pushViewport(vplayout(3:5, 5))\n"
        "}\n"
        "draw_legend(color, breaks, legend, fontsize = fontsize, ...)\n"
        "upViewport()\n"
        "}\n"
        "\n"
        "\n"
        "}\n"
        "\n"
        "generate_breaks = function(x, n){\n"
        "seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)\n"
        "}\n"
        "\n"
        "scale_vec_colours = function(x, col = rainbow(10), breaks = NA){\n"
        "return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])\n"
        "}\n"
        "\n"
        "scale_colours = function(mat, col = rainbow(10), breaks = NA){\n"
        "mat = as.matrix(mat)\n"
        "return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))\n"
        "}\n"
        "\n"
        "cluster_mat = function(mat, distance, method){\n"
        "if(!(method %in% c(\"ward\", \"single\", \"complete\", \"average\", \"mcquitty\", \"median\", \"centroid\"))){\n"
        "stop(\"clustering method has to one form the list: 'ward', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.\")\n"
        "}\n"
        "if(!(distance[1] %in% c(\"correlation\", \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\", \"minkowski\")) & class(distance) != \"dist\"){\n"
        "print(!(distance[1] %in% c(\"correlation\", \"euclidean\", \"maximum\", \"manhattan\", \"canberra\", \"binary\", \"minkowski\")) | class(distance) != \"dist\")\n"
        "stop(\"distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'\")\n"
        "}\n"
        "if(distance[1] == \"correlation\"){\n"
        "d = dist(1 - cor(t(mat)))\n"
        "}\n"
        "else{\n"
        "if(class(distance) == \"dist\"){\n"
        "d = distance\n"
        "}\n"
        "else{\n"
        "d = dist(mat, method = distance)\n"
        "}\n"
        "}\n"
        "\n"
        "return(hclust(d, method = method))\n"
        "}\n"
        "\n"
        "scale_rows = function(x){\n"
        "m = apply(x, 1, mean, na.rm = T)\n"
        "s = apply(x, 1, sd, na.rm = T)\n"
        "return((x - m) / s)\n"
        "}\n"
        "\n"
        "scale_mat = function(mat, scale){\n"
        "if(!(scale %in% c(\"none\", \"row\", \"column\"))){\n"
        "stop(\"scale argument shoud take values: 'none', 'row' or 'column'\")\n"
        "}\n"
        "mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))\n"
        "return(mat)\n"
        "}\n"
        "\n"
        "generate_annotation_colours = function(annotation, annotation_colors, drop){\n"
        "if(is.na(annotation_colors)[[1]][1]){\n"
        "annotation_colors = list()\n"
        "}\n"
        "count = 0\n"
        "for(i in 1:ncol(annotation)){\n"
        "if(is.character(annotation[, i]) | is.factor(annotation[, i])){\n"
        "if (is.factor(annotation[, i]) & !drop){\n"
        "count = count + length(levels(annotation[, i]))\n"
        "}\n"
        "count = count + length(unique(annotation[, i]))\n"
        "}\n"
        "}\n"
        "\n"
        "factor_colors = hsv((seq(0, 1, length.out = count + 1)[-1] +\n"
        "0.2)%%1, 0.7, 0.95)\n"
        "\n"
        "set.seed(3453)\n"
        "\n"
        "for(i in 1:ncol(annotation)){\n"
        "if(!(colnames(annotation)[i] %in% names(annotation_colors))){\n"
        "if(is.character(annotation[, i]) | is.factor(annotation[, i])){\n"
        "n = length(unique(annotation[, i]))\n"
        "if (is.factor(annotation[, i]) & !drop){\n"
        "n = length(levels(annotation[, i]))\n"
        "}\n"
        "ind = sample(1:length(factor_colors), n)\n"
        "annotation_colors[[colnames(annotation)[i]]] = factor_colors[ind]\n"
        "l = levels(as.factor(annotation[, i]))\n"
        "l = l[l %in% unique(annotation[, i])]\n"
        "if (is.factor(annotation[, i]) & !drop){\n"
        "l = levels(annotation[, i])\n"
        "}\n"
        "names(annotation_colors[[colnames(annotation)[i]]]) = l\n"
        "factor_colors = factor_colors[-ind]\n"
        "}\n"
        "else{\n"
        "r = runif(1)\n"
        "annotation_colors[[colnames(annotation)[i]]] = hsv(r, c(0.1, 1), 1)\n"
        "}\n"
        "}\n"
        "}\n"
        "return(annotation_colors)\n"
        "}\n"
        "\n"
        "kmeans_pheatmap = function(mat, k = min(nrow(mat), 150), sd_limit = NA, ...){\n"
        "# Filter data\n"
        "if(!is.na(sd_limit)){\n"
        "s = apply(mat, 1, sd)\n"
        "mat = mat[s > sd_limit, ]\n"
        "}\n"
        "\n"
        "# Cluster data\n"
        "set.seed(1245678)\n"
        "km = kmeans(mat, k, iter.max = 100)\n"
        "mat2 = km$centers\n"
        "\n"
        "# Compose rownames\n"
        "t = table(km$cluster)\n"
        "rownames(mat2) = sprintf(\"cl%s_size_%d\", names(t), t)\n"
        "\n"
        "# Draw heatmap\n"
        "pheatmap(mat2, ...)\n"
        "}\n"
        "\n"
        "pheatmap = function(mat, color = colorRampPalette(rev(c(\"#D73027\", \"#FC8D59\", \"#FEE090\", \"#FFFFBF\", \"#E0F3F8\", \"#91BFDB\", \"#4575B4\")))(100), kmeans_k = NA, breaks = NA, border_color = \"grey60\", cellwidth = NA, cellheight = NA, scale = \"none\", cluster_rows = TRUE, cluster_cols = TRUE, clustering_distance_rows = \"euclidean\", clustering_distance_cols = \"euclidean\", clustering_method = \"complete\",  treeheight_row = ifelse(cluster_rows, 50, 0), treeheight_col = ifelse(cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA, legend_labels = NA, annotation = NA, annotation_colors = NA, annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, number_format = \"%.2f\", fontsize_number = 0.8 * fontsize, filename = NA, width = NA, height = NA, ...){\n"
        "\n"
        "# Preprocess matrix\n"
        "mat = as.matrix(mat)\n"
        "mat = scale_mat(mat, scale)\n"
        "\n"
        "# Kmeans\n"
        "if(!is.na(kmeans_k)){\n"
        "# Cluster data\n"
        "km = kmeans(mat, kmeans_k, iter.max = 100)\n"
        "mat = km$centers\n"
        "\n"
        "# Compose rownames\n"
        "t = table(km$cluster)\n"
        "rownames(mat) = sprintf(\"cl%s_size_%d\", names(t), t)\n"
        "}\n"
        "else{\n"
        "km = NA\n"
        "}\n"
        "\n"
        "# Do clustering\n"
        "if(cluster_rows){\n"
        "tree_row = cluster_mat(mat, distance = clustering_distance_rows, method = clustering_method)\n"
        "mat = mat[tree_row$order, ]\n"
        "}\n"
        "else{\n"
        "tree_row = NA\n"
        "treeheight_row = 0\n"
        "}\n"
        "\n"
        "if(cluster_cols){\n"
        "tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, method = clustering_method)\n"
        "mat = mat[, tree_col$order]\n"
        "}\n"
        "else{\n"
        "tree_col = NA\n"
        "treeheight_col = 0\n"
        "}\n"
        "\n"
        "# Format numbers to be displayed in cells\n"
        "if(display_numbers){\n"
        "fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))\n"
        "attr(fmat, \"draw\") = TRUE\n"
        "}\n"
        "else{\n"
        "fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))\n"
        "attr(fmat, \"draw\") = FALSE\n"
        "}\n"
        "\n"
        "\n"
        "# Colors and scales\n"
        "if(!is.na(legend_breaks[1]) & !is.na(legend_labels[1])){\n"
        "if(length(legend_breaks) != length(legend_labels)){\n"
        "stop(\"Lengths of legend_breaks and legend_labels must be the same\")\n"
        "}\n"
        "}\n"
        "\n"
        "\n"
        "if(is.na(breaks[1])){\n"
        "breaks = generate_breaks(as.vector(mat), length(color))\n"
        "}\n"
        "if (legend & is.na(legend_breaks[1])) {\n"
        "legend = grid.pretty(range(as.vector(breaks)))\n"
        "names(legend) = legend\n"
        "}\n"
        "else if(legend & !is.na(legend_breaks[1])){\n"
        "legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]\n"
        "\n"
        "if(!is.na(legend_labels[1])){\n"
        "legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]\n"
        "names(legend) = legend_labels\n"
        "}\n"
        "else{\n"
        "names(legend) = legend\n"
        "}\n"
        "}\n"
        "else {\n"
        "legend = NA\n"
        "}\n"
        "mat = scale_colours(mat, col = color, breaks = breaks)\n"
        "\n"
        "# Preparing annotation colors\n"
        "if(!is.na(annotation[[1]][1])){\n"
        "annotation = annotation[colnames(mat), , drop = F]\n"
        "annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)\n"
        "}\n"
        "\n"
        "if(!show_rownames){\n"
        "rownames(mat) = NULL\n"
        "}\n"
        "\n"
        "if(!show_colnames){\n"
        "colnames(mat) = NULL\n"
        "}\n"
        "\n"
        "# Draw heatmap\n"
        "heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, ...)\n"
        "garbage<-dev.off()\n"
        "invisible(list(tree_row = tree_row, tree_col = tree_col, kmeans = km))\n"
        "\n"
        "}\n"
        "\n"
        "\n"
        "\n"
        "\n"
        "\n"
        "\n"
        "\n"
        "filename=\""<<matrix_filename_str<<"\"\n"
        "d <- read.table(filename, header=TRUE, sep=\"\\t\", check.names=FALSE)\n"
        "#d <- d[order(d$.),]\n"
        "row.names(d) <- d$.\n"
        "d <- d[,3:ncol(d)-1]\n"
        "dmat<-data.matrix(d)\n"
        "\n"
        "mycol=colorRampPalette(c(\"white\",\"yellow\",\"#FF9B0D\",\"#BD2400\",\"#631300\",\"black\"))\n"
        "maintitle=tail(strsplit(strsplit(filename, '.fastq')[[1]][1],'/')[[1]],1)\n"
        "maintitle=tail(strsplit(strsplit(maintitle, '.fastq')[[1]][1],'\\\\\\\\')[[1]],1)\n"

        "pheatmap(dmat, cellheight=10, cellwidth=10, cluster_rows=FALSE, filename=paste(filename, '.pdf', sep=''), cluster_cols=FALSE, legend=TRUE, main=maintitle, color=mycol(75), breaks=c(seq(0,0.75,0.01)))\n"
        "\n";

        R.close();
        std::string Rcommand = "Rscript "+R_filename_str;
        int result=::system(Rcommand.c_str());
        remove(R_file);
        if(result==-1) {
            std::cout<<"Error while launching Rscript"<<std::endl;
            return -1;
        }
        else{
            return 0;
        }


}


int R_quality(char* directory, std::string filepath, std::string filename, std::string quality_filename_str, int extra_columns){

    const char* R_file;
    std::string R_filename_str;
    if(directory){
        R_filename_str=directory+filename+".R";
        R_file=R_filename_str.c_str();
    }else{
        R_filename_str=filepath+".R";
        R_file=R_filename_str.c_str();
    }

    std::ofstream R (R_file, std::ofstream::out);
    if(not R.is_open()){
        std::cout<<"Error: Failure opening "<<R_file<<" for writing."<<std::endl;
        return 1;
    }
    for (int i = 0; i < quality_filename_str.length(); ++i)
        if (quality_filename_str[i] == '\\')
        {
        quality_filename_str.insert(i, 1, '\\');
        ++i; // Skip inserted char
        }
    R<<"library(package=grid)\n"
        "\n"
        "filename=\""<<quality_filename_str<<"\"\n"
        "extra_columns="<<extra_columns<<"\n"
        "d <- read.table(filename, header=TRUE)\n"
        "\n"
        "pdf( paste(filename,'.pdf', sep = '') )\n"
        "par(mar=c(5,6,6,2) + 0.1, oma=c(3,0,0,0), mgp = c(4, 1, 0))\n"
        "\n"
        "varx=1:length(d[,1])\n"
        "vary=d[,1]\n"
        "plot(x=varx, y=vary, col='red', ylim=c(0,0.631), type='l', lwd=3, las=1, xlab='Position along read',col.lab=rgb(0,0.5,0), ylab='Mean probability of error', col.lab=rgb(0,0.5,0))\n"
        "abline(h=seq(0.1,0.6,0.1), lty=1, col='grey')\n"
        "\n"
        "legend(0.59,c(\"Error across all tiles\", \"Individual base call error per tile\"), cex=0.7, col=c('red', 'black'), lty=c(1,3), lwd=c(3,1))\n"
        "Title=tail(strsplit(strsplit(filename, '.fastq')[[1]][1],'/')[[1]],1)\n"
        "Title=tail(strsplit(strsplit(Title, '.fastq')[[1]][1],'\\\\\\\\')[[1]],1)\n"
        "title(paste(\'Sample: \', Title))\n"
        "mtext(\"Average base call error per nucleotide position\", 3, line=1)\n"
        "\n"
        "for(i in c(1:length(d[1,]))){\n"
        "\n"
        "if( ( i + extra_columns ) %% ( 1 + extra_columns ) == 0 ){\n"
        "\n"
        "lines(1:length(d[,i]), d[,i], lty=3, col='black')\n"
        "}\n"
        "}\n"
        "\n"
        "lines(1:length(d[,1]), d[,1], type='l', col='red')\n"
        "garbage<-dev.off()\n";
    R.close();
    std::string Rcommand = "Rscript "+R_filename_str;
    int result=::system(Rcommand.c_str());
    remove(R_file);
    if(result==-1) {
        std::cout<<"Error while launching Rscript"<<std::endl;
        return -1;
    }
    else{
        return 0;
    }

}


int R_histogram(char* directory, std::string filepath, std::string filename, std::string segments_filename_str, double prob_cutoff){

    const char* R_file;
    std::string R_filename_str;
    if(directory){
        R_filename_str=directory+filename+".R";
        R_file=R_filename_str.c_str();
    }else{
        R_filename_str=filepath+".R";
        R_file=R_filename_str.c_str();
    }

    std::ofstream R (R_file, std::ofstream::out);
    if(not R.is_open()){
        std::cout<<"Error: Failure opening "<<R_file<<" for writing."<<std::endl;
        return 1;
    }
    for (int i = 0; i < segments_filename_str.length(); ++i)
        if (segments_filename_str[i] == '\\')
        {
        segments_filename_str.insert(i, 1, '\\');
        ++i; // Skip inserted char
        }
    R<<"filename <- \""<<segments_filename_str<<"\"\n"
        "cutoff <- "<<prob_cutoff<<"\n"
        "output <- \""<<segments_filename_str<<"_hist\"\n"
        "\n"
        "d <- read.table(filename, header=T)\n"
        "maxup=max(d[,2])\n"
        "pdf( paste(output,'.pdf', sep = ''), width = 11 )\n"
        "if (maxup>0.45){\n"
        "# I want to plot the lower values up to 55, then a split to 95 for the\n"
        "# last top. This should make it clear which is the highest, without\n"
        "# drowning out the other data.\n"
        "\n"
        "# I want the split to be approx 5% of the scale,\n"
        "\n"
        "# as I am to plot the ranges 0 - 55 and 95 - 140, in total 10 decades,\n"
        "lower=c(0,0.4)\n"
        "upper=c(maxup-0.05, maxup)\n"
        "# This is 10 decades. I multiply that with 2 and add 5% and get 21 units on the outer\n"
        "# Y axis:\n"
        "y_outer=(lower[2]+upper[1]-upper[2])*100\n"
        "lowspan=c(0,(2*y_outer/3)-1)\n"
        "topspan=c(lowspan[2]+5, y_outer)\n"
        "\n"
        "\n"
        "cnvrt.coords <-function(x,y=NULL){\n"
        "# Stolen from the teachingDemos library, simplified for this use case\n"
        "xy <- xy.coords(x,y, recycle=TRUE)\n"
        "cusr <- par('usr')\n"
        "cplt <- par('plt')\n"
        "plt <- list()\n"
        "plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])\n"
        "plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])\n"
        "fig <- list()\n"
        "fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]\n"
        "fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]\n"
        "return( list(fig=fig) )\n"
        "}\n"
        "\n"
        "subplot <- function(fun, x, y=NULL){\n"
        "# Stolen from the teachingDemos l	ibrary, simplified for this use case\n"
        "old.par <- par(no.readonly=TRUE)\n"
        "on.exit(par(old.par))\n"
        "xy <- xy.coords(x,y)\n"
        "xy <- cnvrt.coords(xy)$fig\n"
        "par(plt=c(xy$x,xy$y), new=TRUE)\n"
        "fun\n"
        "tmp.par <- par(no.readonly=TRUE)\n"
        "return(invisible(tmp.par))\n"
        "}\n"
        "\n"
        "##############################################\n"
        "#\n"
        "#\n"
        "# The main program starts here:\n"
        "#\n"
        "#\n"
        "\n"
        "# Setting up an outer wireframe for the plots.\n"
        "par(mar=c(8,6,6,3) + 0.1, oma=c(0,0,0,0), mgp = c(3, 1, 0))\n"
        "plot(c(0,1),c(0,y_outer),type='n',axes=FALSE,xlab=paste(\'Length of longest contiguous read segments with quality higher than\',cutoff),col.lab=rgb(0,0.5,0), ylab='Proportion of reads', col.lab=rgb(0,0.5,0))\n"
        "Title=tail(strsplit(strsplit(filename, '.fastq')[[1]][1],'/')[[1]],1)\n"
        "Title=tail(strsplit(strsplit(Title, '.fastq')[[1]][1],'\\\\\\\\')[[1]],1)\n"
        "title(paste(\'Sample: \',Title))\n"
        "mtext(paste(\"p cutoff = \",cutoff), 3, line=1)\n"
        "mtext(\"Sum of the segments = 1\", 1, line=6)\n"
        "# Plotting the lower range in the lower 11/21 of the plot.\n"
        "# xpd=FALSE to clip the bars\n"
        "tmp<-subplot(barplot(d[,2],col=\'blue\',ylim=lower,xpd=FALSE,las=0, names.arg = d[,1], las=1), x=c(0,1),y=lowspan)\n"
        "op <- par(no.readonly=TRUE)\n"
        "par(tmp)\n"
        "abline(h=seq(0.0,0.4,0.1), lty=1, col=\'grey\')\n"
        "par(op)\n"
        "subplot(barplot(d[,2],col=\'blue\',ylim=lower,xpd=FALSE,las=0, names.arg = d[,1], las=1), x=c(0,1),y=lowspan)\n"
        "\n"
        "\n"
        "# Plotting the upper range in the upper 9/21 of the plot, 1/21 left to\n"
        "# the split. Again xpd=FALSE, names.arg is set up to avoid having\n"
        "# the names plotted here, must be some easier way to do this but\n"
        "# this works\n"
        "\n"
        "tmp<-subplot(barplot(d[,2],col=\'blue\',ylim=c(round(upper[1], digits = 1),round(upper[2], digits = 1)+0.1), xpd=FALSE, las=1), x=c(0,1),y=topspan)\n"
        "op <- par(no.readonly=TRUE)\n"
        "par(tmp)\n"
        "abline(h=seq(round(upper[1], digits = 1),round(upper[2], digits = 1)+0.1,0.1), lty=1, col=\'grey\')\n"
        "par(op)\n"
        "subplot(barplot(d[,2],col=\'blue\',ylim=c(round(upper[1], digits = 1),round(upper[2], digits = 1)+0.1), xpd=FALSE, las=1), x=c(0,1),y=topspan)\n"
        "\n"
        "# Legend. An annoiance is that the colors comes in the opposite\n"
        "# order than in the plot.\n"
        "\n"
        "#legend(0.05,26,c(\'Reads\'), cex=0.7, col=c(\'blue\'), pch=15)\n"
        "\n"
        "# so far so good. (Just run the upper part to see the result so far)\n"
        "# Just want to make the ends of the axes a bit nicer.\n"
        "# All the following plots are in units of the outer coordinate system\n"
        "\n"
        "lowertop=lowspan[2]+(topspan[1]-lowspan[2])/2  # Where to end the lower axis\n"
        "breakheight=1   # Height of the break\n"
        "upperbot=lowertop+breakheight#(lowspan[2]+(topspan[1]-lowspan[2])/2)+breakheight# Where to start the upper axes\n"
        "markerheight=0.5 # Heightdifference for the break markers\n"
        "markerwidth=.03  # With of the break markers\n"
        "\n"
        "# Draw the break markers:\n"
        "#lines(c(0,0),c(1,lowertop))\n"
        "lines(c(markerwidth/-2,markerwidth/2),c(lowertop-markerheight/2,lowertop+markerheight/2))\n"
        "#lines(c(0,0),c(upperbot-breakheight,14))\n"
        "#lines(c(0,0),c(upperbot,maxup))\n"
        "lines(c(markerwidth/-2,markerwidth/2),c(upperbot-markerheight/2,upperbot+markerheight/2))\n"
        "\n"
        "}else{\n"
        "\n"
        "par(mar=c(8,6,6,3) + 0.1, oma=c(0,0,0,0), mgp = c(3, 1, 0))\n"
        "barplot(d[,2], names.arg = d[,1], space = 0, ylim=c(0, 0.4), col=\'blue\', las=1, xlab=paste(\'Length of longest contiguous read segments with quality higher than\',cutoff),col.lab=rgb(0,0.5,0), ylab=\'Proportion of reads\', col.lab=rgb(0,0.5,0), axis.lty = 1, cex.names = 0.9 )\n"
        "abline(h=seq(0.1,0.4,0.1), lty=1, col=\'grey\')\n"
        "barplot(d[,2], add=TRUE, names.arg = d[,1], space = 0, ylim=c(0, 0.4), col=\'blue\', las=1, xlab=paste(\'Length of longest contiguous read segments with quality higher than\',cutoff),col.lab=rgb(0,0.5,0), ylab=\'Proportion of reads\', col.lab=rgb(0,0.5,0), axis.lty = 1, cex.names = 0.9 )\n"
        "\n"
        "Title=tail(strsplit(strsplit(filename, '.fastq')[[1]][1],'/')[[1]],1)\n"
        "Title=tail(strsplit(strsplit(Title, '.fastq')[[1]][1],'\\\\\\\\')[[1]],1)\n"
        "title(paste(\'Sample: \',Title))\n"
        "mtext(paste(\"p cutoff = \",cutoff), 3, line=1)\n"
        "mtext(\"Sum of the segments = 1\", 1, line=6)\n"
        "\n"
        "#legend(0.05,0.35,c(\'Reads\'), cex=0.7, col=c(\'blue\'), pch=15)\n"
        "}\n"
        "garbage<-dev.off()\n";

    R.close();
    std::string Rcommand = "Rscript "+R_filename_str;
    int result=::system(Rcommand.c_str());
    remove(R_file);
    if(result==-1) {
        std::cout<<"Error while launching Rscript"<<std::endl;
        return -1;
    }
    else{
        return 0;
    }

}


int R_cumulative(char* directory, std::string filepath, std::string filename, std::string segments_filename_str, double prob_cutoff){
    const char* R_file;
    std::string R_filename_str;

    for (int i = 0; i < segments_filename_str.length(); ++i)
    if (segments_filename_str[i] == '\\')
    {
    segments_filename_str.insert(i, 1, '\\');
    ++i; // Skip inserted char
    }

    if(directory){
        R_filename_str=directory+filename+".R";
        R_file=R_filename_str.c_str();
    }else{
        R_filename_str=filepath+".R";
        R_file=R_filename_str.c_str();
    }

    std::ofstream R (R_file, std::ofstream::out);
    if(not R.is_open()){
        std::cout<<"Error: Failure opening "<<R_file<<" for writing."<<std::endl;
        return 1;
    }
    R<<"filename <- \""<<segments_filename_str<<"\"\n"
        "output <- \""<<segments_filename_str<<"_cumulative\"\n"
        "cutoff <- "<<prob_cutoff<<"\n"
        "\n"
        "d<-read.table(filename, header=TRUE)\n"
        "d$sumcol=d$proportion_of_reads\n"
        "\n"
        "#d$sumcol=d$sumcol+d$proportion_of_reads\n"
        "sum=0\n"
        "\n"
        "d$ideal=1\n"
        "d$ideal[nrow(d)]=0\n"
        "for(i in 1:nrow(d)){\n"
        "sum=sum+d$proportion_of_reads[i]\n"
        "d$sumcol[i]=1-sum}\n"
        "\n"
        "pdf( paste(output,'.pdf', sep = ''), width = 11 )\n"
        "par(mar=c(4.5,6,6,2) + 0.1, oma=c(0,0,0,0), mgp = c(3, 1, 0))\n"
        "plot(x=d$read_length, y=d$sumcol, panel.first = abline(h=seq(0.0,1,0.2), lty=1, col='grey'),type='l', xlab=paste(\'Length of longest contiguous read segments with quality higher than\',cutoff),col.lab=rgb(0,0.5,0), ylab='Proportion of reads',col.lab=rgb(0,0.5,0), lwd=2, col='purple', las=1)\n"
        "lines(d$read_length,d$ideal, type='l', lty=3, col='black', lwd=2)\n"
        "\n"
        "\n"
        "legend(0.77,c(\"Reads\", \"Ideal\"), lty=c(1,3), cex=0.7, col=c('purple', 'black'), lwd=2)\n"
        "Title=tail(strsplit(strsplit(filename, '.fastq')[[1]][1],'/')[[1]],1)\n"
        "Title=tail(strsplit(strsplit(Title, '.fastq')[[1]][1],'\\\\\\\\')[[1]],1)\n"
        "title(paste(\'Sample: \',Title))\n"
        "mtext(paste('p cutoff = ',cutoff), 3, line=1)\n"
        "garbage<-dev.off()\n";
    R.close();
    std::string Rcommand = "Rscript "+R_filename_str;
    int result=::system(Rcommand.c_str());
    remove(R_file);
    remove("Rplots.pdf");

    if(result==-1) {
        std::cout<<"Error while launching Rscript"<<std::endl;
        return -1;
    }
    else{
        return 0;
    }

}


int R_length_summary(char* directory, std::string filepath, std::string filename, std::string summary_filename_str, bool paired, int length){
    const char* R_file;
    std::string R_filename_str;
    if(directory){
        R_filename_str=directory+filename+".R";
        R_file=R_filename_str.c_str();
    }else{
        R_filename_str=filepath+".R";
        R_file=R_filename_str.c_str();
    }
    for (int i = 0; i < summary_filename_str.length(); ++i)
        if (summary_filename_str[i] == '\\')
        {
        summary_filename_str.insert(i, 1, '\\');
        ++i; // Skip inserted char
        }
    std::ofstream R (R_file, std::ofstream::out);
    if(not R.is_open()){
        std::cout<<"Error: Failure opening "<<R_file<<" for writing."<<std::endl;
        return 1;
    }
    if(paired==true){
        R<<"options(scipen=5)\n"
            "length="<<length<<"\n"
            "filename <- \""<<summary_filename_str<<"\"\n"
            "output <- \""<<summary_filename_str<<"\"\n"
            "d<-read.table(filename, header=FALSE)\n"
            "dnew<-cbind(c(\"Paired\n(Read1 and Read2)\", \"Single\", \"Discard\"), c(d$V2[1]+d$V2[2], d$V2[3], d$V2[4]))\n"
            "uplimit<-max(as.integer(dnew[,2]))\n"
            "total=sum(d$V2)\n"

            "pdf( paste(output, '.pdf', sep = ''))\n"
            "par(mar=c(4.5,6,6,2) + 0.1, oma=c(0,0,0,0), mgp = c(3, 1, 0), cex.lab=1, cex.main=1.3)\n"
            "barX<-barplot(as.integer(dnew[,2]), cex.names=0.8, names.arg=(dnew[,1]), font.lab=2, ylab=\"Reads\", ylim=c(0,uplimit+uplimit/100*20), col=c(\"darkgreen\", \"dodgerblue4\", \"firebrick3\"), beside=TRUE, main=\"Sorted reads\")\n"
            "mtext(paste('sorted length = ',length), 3, line=1)\n"
            "text(cex=0.9, x=barX, y=as.integer(dnew[,2])+par(\"cxy\")[2]/2, paste(round(as.integer(dnew[,2]),2), \" (\", sprintf(\"%.2f\", as.integer(dnew[,2])*100/total), \"%)\", sep =''), xpd=TRUE, col=c(\"darkgreen\", \"dodgerblue4\", \"firebrick3\"))\n"
            "garbage<-dev.off()\n";
    }else{
        R<<"options(scipen=5)\n"
            "length="<<length<<"\n"
            "filename <- \""<<summary_filename_str<<"\"\n"
            "output <- \""<<summary_filename_str<<"\"\n"
            "d<-read.table(filename, header=FALSE)\n"
            "dnew<-cbind(c(\"Single\", \"Discard\"), c(d$V2[1], d$V2[2]))\n"
            "uplimit<-max(as.integer(dnew[,2]))\n"
            "total=sum(d$V2)\n"

            "pdf( paste(output, '.pdf', sep = ''))\n"
            "par(mar=c(4.5,6,6,2) + 0.1, oma=c(0,0,0,0), mgp = c(3, 1, 0), cex.lab=1, cex.main=1.3)\n"
            "barX<-barplot(as.integer(dnew[,2]), names.arg=(dnew[,1]), font.lab=2, ylab=\"Reads\", ylim=c(0,uplimit+uplimit/100*20), col=c(\"dodgerblue4\", \"firebrick3\"), beside=TRUE, main=\"Sorted reads\")\n"
            "mtext(paste('sorted length = ',length), 3, line=1)\n"
            "text(cex=0.9, x=barX, y=as.integer(dnew[,2])+par(\"cxy\")[2]/2, paste(round(as.integer(dnew[,2]),2), \" (\", sprintf(\"%.2f\", as.integer(dnew[,2])*100/total), \"%)\", sep =''), xpd=TRUE, col=c(\"dodgerblue4\", \"firebrick3\"))\n"
            "garbage<-dev.off()\n";
    }

    R.close();
    std::string Rcommand = "Rscript "+R_filename_str;
    int result=::system(Rcommand.c_str());
    remove(R_file);
    remove("Rplots.pdf");

    if(result==-1) {
        std::cout<<"Error while launching Rscript"<<std::endl;
        return -1;
    }
    else{
        return 0;
    }


}
#endif
