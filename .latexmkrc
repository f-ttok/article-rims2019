#!/usr/bin/env perl
$pdf_mode = 3;
$latex = "uplatex -interaction=nonstopmode -synctex=1";
$bibtex = "pbibtex";
$dvipdf = "dvipdfmx %O -o %D %S";
$out_dir = "./output";
$max_repeat = 2;
$pvc_view_file_via_temporary = 0;
