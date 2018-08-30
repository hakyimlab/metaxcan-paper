#! /usr/bin/env python
import os
import subprocess
import svgwrite
import math
import shutil

########################################################################################################################
def ensure_requisite_folders(path):
    folder = os.path.split(path)[0]
    if len(folder) and not os.path.exists(folder):
        os.makedirs(folder)

def _png_name(p):
    return p.split(".svg")[0]+".png"

def to_png(from_path, to_path):
    ensure_requisite_folders(to_path)
    cmd = \
"""
convert {} {}
""".format(from_path, to_path)
    subprocess.call(cmd.split())

def _advance_cursor(c, x, y):
    return (c[0]+x, c[1]+y)

def _label(dwg, _c, text, offset=(0,0),
           style="font-size:40;font-family:Arial;font-weight:bold;stroke:black;stroke-width:1;fill:black"):
    dwg.add(dwg.text(text, insert=(_c[0] + offset[0], _c[1]+offset[1]), fill='black',style=style))

########################################################################################################################
def fig_1(args):
    _p = "_fig1.svg"#os.path.join(args.output_folder, "fig1.svg")

    #Rationale:
    # Fig1_predixcan_gwas_summary_predixcan.png is 891x639, population comparison is 1200 * 1200, WT is 600 * 1200
    _1_size = (600, 600)
    _2_size = (300, 600)
    _size = (50+_1_size[0]+50+_1_size[0]+50+_2_size[0], _1_size[1])

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    _c = (50, 0) # conceptual cursor
    dwg.add(dwg.image(os.path.join(args.input_folder, "predixcan_metaxcan_simulated_comparison.png"), _c, _1_size))
    _label(dwg, _c, "a", (-25,50))

    _c = _advance_cursor(_c, _1_size[1]+50,0)
    dwg.add(dwg.image(os.path.join(args.input_folder,"predixcan_metaxcan_igrowth_comparison.png"), _c, _1_size))
    _label(dwg, _c, "b", (-25,50))

    _c =_advance_cursor (_c, _1_size[0]+50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder,"predixcan_metaxcan_t1d_bd_comparison.png"), _c, _2_size))
    _label(dwg, _c, "c", (-25, 50))

    dwg.save()
    t = os.path.join(args.output_folder, "fig-s-predixcan-performance.png")
    to_png(_p, t)
    os.remove(_p)

#
def fig_2(args):
    shutil.copy(os.path.join(args.input_folder, "clinvar.png"),
                os.path.join(args.output_folder, "fig-clinvar.png"))

    shutil.copy(os.path.join(args.input_folder, "clinvar_hla_excluded.png"),
                os.path.join(args.output_folder, "supp-fig-clinvar-hla-excluded.png"))

    #shutil.copy(os.path.join(args.input_folder_2, "metaxcan_application.png"),
    #            os.path.join(args.output_folder, "fig-metaxcan-application.png"))

    #shutil.copy(os.path.join(args.input_folder_2, "metaxcan_framework.png"),
    #            os.path.join(args.output_folder, "fig-metaxcan-framework.png"))

    shutil.copy(os.path.join(args.input_folder, "replication_ldl.png"),
                os.path.join(args.output_folder, "fig-replication-ldl.png"))

    shutil.copy(os.path.join(args.input_folder, "replication_nodiagram.png"),
                os.path.join(args.output_folder, "supp-fig-replication-gera.png"))

    shutil.copy(os.path.join(args.input_folder, "geuvadis_expression_e_short.png"),
                os.path.join(args.output_folder, "supp-fig-geuvadis-expression-short.png"))

    shutil.copy(os.path.join(args.input_folder, "predixcan_metaxcan_alt_comparison.png"),
               os.path.join(args.output_folder, "supp-fig-predixcan-metaxcan-alt-comparison.png"))
               
    shutil.copy(os.path.join(args.input_folder, "twas_predixcan_nsignificant.png"),
                os.path.join(args.output_folder, "supp-fig-twas-predixcan-nsignificant.png"))
    # Yeah, no.
    # shutil.copy(os.path.join(args.input_folder, "predixcan_metaxcan_simulated_alt_comparison.png"),
    #            os.path.join(args.output_folder, "supp-fig-predixcan-metaxcan-simulated_alt_comparison.png"))

    shutil.copy(os.path.join(args.input_folder, "table-absZ2-tissue-phenotype-no-high.png"),
                os.path.join(args.output_folder, "supp-fig-association-enrichment.png"))

#
def fig_3(args):
    _p = "_fig3.svg"
    #twas image is 1404*1524, twas vs spredixcan is 500, twas vs spredixcan colocalization is 700
    _1_size = (math.ceil(500.0/1524*1404), 500)
    _2_size = (500, 500)
    _size = (50+_2_size[0]+50+_2_size[0], _1_size[1]+50+_2_size[1])

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    _c = (50+ math.ceil((_2_size[0]-_1_size[0])/2), 0) # conceptual cursor
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "s_predixcan_twas.png"), _c, _1_size))
    _label(dwg, (25,50), "a")

    _c = (_size[0] -_2_size[0], 0)
    dwg.add(dwg.image(os.path.join(args.input_folder,"s_predixcan_vs_twas.png"), _c, _2_size))
    _label(dwg, _c, "b", (-25,50))

    _c =(50, _1_size[1]+50)
    dwg.add(dwg.image(os.path.join(args.input_folder,"twas_predixcan_p3_proportion.png"), _c, _2_size))
    _label(dwg, _c, "c", (-25,50))

    _c = _advance_cursor(_c, _2_size[0]+50,0)
    dwg.add(dwg.image(os.path.join(args.input_folder,"twas_predixcan_p4_proportion.png"), _c, _2_size))
    _label(dwg, _c, "d", (-25, 50))

    dwg.save()
    t = os.path.join(args.output_folder, "fig-twas-comparison.png")
    to_png(_p, t)
    os.remove(_p)

#
def fig_4_deprecated(args):
    _p = "_fig4.svg"
    #smr image is 703x1093
    _1_size = (500.0/1093*703, 500)
    _2_size = (500, 500)
    _size = ((50 + 500)*3, 500 + 50 + 500)

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    _c = ((550-_1_size[0])/2,0) #conceptual cursor
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "smr.png"), _c, _1_size))
    _label(dwg, _c, "a", (-25,50))

    _c = (600, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "smr_simul_2.png"), _c, _2_size))
    _label(dwg, _c, "b", (-25, 50))

    _c = _advance_cursor(_c, _2_size[0]+50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "smr_meta_simul_2.png"), _c, _2_size))
    _label(dwg, _c, "c", (-25, 50))

    _c = (50, 550)
    dwg.add(dwg.image(os.path.join(args.input_folder, "s_predixcan_vs_smr.png"), _c, _2_size ))
    _label(dwg, _c, "d", (-25, 50))

    _c = _advance_cursor(_c, _2_size[0] + 50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "giant_height_adipose_subcutaneous_smr_vs_eqtl.png"), _c, _2_size))
    _label(dwg, _c, "e", (-25, 50))

    _c = _advance_cursor(_c, _2_size[0] + 50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "giant_height_adipose_subcutaneous_smr_vs_gwas.png"), _c, _2_size))
    _label(dwg, _c, "f", (-25, 50))

    dwg.save()
    t = os.path.join(args.output_folder, "fig-smr-comparison.png")
    to_png(_p, t)
    os.remove(_p)

def fig_4(args):
    _p = "_fig4.svg"
    #smr image is 703x1093
    _1_size = (500.0/1093*703, 500)
    _2_size = (500, 500)
    _size = ((50 + 500)*4, 500 + 50 + 500)

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    _c = ((550-_1_size[0])/2,0) #conceptual cursor
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "smr.png"), _c, _1_size))
    _label(dwg, _c, "a", (-25,50))

    _c = (600, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "s_predixcan_vs_smr.png"), _c, _2_size))
    _label(dwg, _c, "b", (-25,50))

    _c = _advance_cursor(_c, _2_size[0]+50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "giant_height_adipose_subcutaneous_smr_vs_eqtl.png"), _c, _2_size))
    _label(dwg, _c, "c", (-25,50))

    _c = _advance_cursor(_c, _2_size[0] + 50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "giant_height_adipose_subcutaneous_smr_vs_gwas.png"), _c, _2_size))
    _label(dwg, _c, "d", (-25, 50))

    _c = (50, 550)
    dwg.add(dwg.image(os.path.join(args.input_folder, "smr_simul_2.png"), _c, _2_size))
    _label(dwg, _c, "e", (-25, 50))

    _c = _advance_cursor(_c, _2_size[0] + 50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "smr_meta_simul_2.png"), _c, _2_size))
    _label(dwg, _c, "f", (-25, 50))

    _c = _advance_cursor(_c, _2_size[0] + 50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "smr_predixcan_p3_proportion.png"), _c, _2_size))
    _label(dwg, _c, "g", (-25, 50))

    _c = _advance_cursor(_c, _2_size[0] + 50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "smr_predixcan_p4_proportion.png"), _c, _2_size))
    _label(dwg, _c, "h", (-25, 50))

    dwg.save()
    t = os.path.join(args.output_folder, "fig-smr-comparison.png")
    to_png(_p, t)
    os.remove(_p)


#
def fig_5(args):
    _p = "_fig5.svg"
    #bubbleplot is 1600*659
    _1_size = (1600, 659)
    _2_size = (1600, 600)
    _size = (_1_size[0]+50, _1_size[1]+50+_2_size[1])

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    _c = (50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "functional_genes_bubbleplot.png"), _c, _1_size))
    _label(dwg, _c, "a", (-25, 50))

    _c = (50, _1_size[1]+50)
    dwg.add(dwg.image(os.path.join(args.input_folder, "significant_results_tissue_specificity_height.png"), _c, _2_size))
    _label(dwg, _c, "b", (-25, 50))

    dwg.save()
    t = os.path.join(args.output_folder, "fig-tissue-specificity.png")
    to_png(_p, t)
    os.remove(_p)

#
def fig_6(args):
    _p = "_fig6.svg"
    # Rationale:
    # fig_colloc_hypothesis.png is 2060 * 2137,
    # gha_eqtl_all_tern_5.png is 800 * 800

    _1_size = (math.ceil(2060.0/2137*700), 700)
    _2_size = (400, 400)
    _size = (_1_size[0]+50+_2_size[0]*2+25, _2_size[1]*2+10)
    _half_y = (_size[1]-_1_size[1])/2

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    _c = (0, _half_y)
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "coloc_hypothesis.png"), _c, _1_size))
    _label(dwg, _c, "a", (10,40))

    _c = (_1_size[0] + 50, 0)
    #dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_p1.png"), _c, _2_size))
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_all_results.png"), _c, _2_size))
    _label(dwg, (_c[0]-10, _half_y+40), "b")

    _c = _advance_cursor(_c, _2_size[0] + 25, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_s.png"), _c, _2_size))
    _label(dwg, (_c[0]-10, _c[1] + _half_y + 40), "c")

    _c = (_1_size[0] + 50, _2_size[1]+10)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_s_heidim05.png"), _c, _2_size))
    _label(dwg, (_c[0]-10, _c[1] + _half_y + 40), "d")

    _c = _advance_cursor(_c, _2_size[0] + 25, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_s_heidiM37.png"), _c, _2_size))
    _label(dwg, (_c[0]-10, _c[1] + _half_y + 40), "e")

    dwg.save()
    t = os.path.join(args.output_folder, "fig-coloc-comparison.png")
    to_png(_p, t)
    os.remove(_p)

#
def fig_supp_coloc(args):
    _p = "_sup_fig_coloc.svg"
    # Rationale:
    # coloc figures are 400*400

    _1_size = (math.ceil(400), 400)
    _size = (_1_size[0]*4+75, _1_size[1]*2+25)
    _half_y = (_size[1]-_1_size[1])/2

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    #first row
    _c = (0, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_all_results.png"), _c, _1_size))
    _label(dwg, _c, "a", (20, 120))

    _c = _advance_cursor(_c, _1_size[0] + 25, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_mpe4.png"), _c, _1_size))
    _label(dwg, _c, "b", (20, 120))

    _c = _advance_cursor(_c, _1_size[0] + 25, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_s.png"), _c, _1_size))
    _label(dwg, _c, "c", (20, 120))

    _c = _advance_cursor(_c, _1_size[0] + 25, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_s_ppq01.png"), _c, _1_size))
    _label(dwg, _c, "d", (20, 120))

    # second row
    _c = (0, _1_size[1] + 25)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_s_heidina.png"), _c, _1_size))
    _label(dwg, _c, "e", (20, 120))

    _c = _advance_cursor(_c, _1_size[0] + 25, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_s_heidim05.png"), _c, _1_size))
    _label(dwg, _c, "f", (20, 120))

    _c = _advance_cursor(_c, _1_size[0] + 25, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_s_heidiM05m37.png"), _c, _1_size))
    _label(dwg, _c, "g", (20, 120))

    _c = _advance_cursor(_c, _1_size[0] + 25, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder, "gha_eqtl_tern_s_heidiM37.png"), _c, _1_size))
    _label(dwg, _c, "h", (20, 120))

    dwg.save()
    t = os.path.join(args.output_folder, "supp-fig-coloc.png")
    to_png(_p, t)
    os.remove(_p)

#
def fig_7(args):
    _p = "_fig7.svg"
    #S-PrediXcan Formula: 959x639; S-PrediXcan to GWAS: 891x639
    _1_size = (891, 639)
    _2_size = (math.ceil(1100.0/850*639), 639)
    _size = (20+_1_size[0]+50+_2_size[0], _2_size[1])

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    _c = (20, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "Fig1_predixcan_gwas_summary_predixcan.png"), _c, _1_size))
    _label(dwg, _c, "a", (0, 40))

    _c = _advance_cursor(_c, _1_size[0] + 50, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "fig-s-predixcan-formula.png"), _c, _1_size))
    _label(dwg, _c, "b", (30, 40))

    dwg.save()
    t = os.path.join(args.output_folder, "fig-s-predixcan.png")
    to_png(_p, t)
    os.remove(_p)

#
def fig_7_vertical(args):
    _p = "_fig7.svg"
    #S-PrediXcan Formula: 1100x850; S-PrediXcan to GWAS: 891x639
    _1_size = (891, 639)
    _2_size = (math.ceil(1100.0 / 850 * 639), 639)
    _size = (_1_size[0]+20, _1_size[1]+10+_2_size[1])

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    _c = (20,0)
    _x = _c[0]
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "Fig1_predixcan_gwas_summary_predixcan.png"), _c, _1_size))
    _label(dwg, (_x,_c[1]+40), "a")

    _c = (math.ceil((_size[0]-_2_size[0])/2.0)+20, _1_size[1]+10)
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "fig-s-predixcan-formula.png"), _c, _2_size))
    _label(dwg, (_x, _c[1]+90), "b")

    dwg.save()
    t = os.path.join(args.output_folder, "fig-s-predixcan.png")
    to_png(_p, t)
    os.remove(_p)

def fig_8_deprecated(args):
  _p = "_fig8.svg"
  _1_size = (700, 700)
  _size = (2*_1_size[0] + 60, _1_size[1])
  
  dwg = svgwrite.Drawing(_p, size=_size)
  dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))
  
  #first row
  _c = (10, 0)
  dwg.add(dwg.image(os.path.join(args.input_folder, "smr_predixcan_p3_proportion.png"), _c, _1_size))
  _label(dwg, _c, "a", (0, 50))

  _c = _advance_cursor(_c, _1_size[0]+50, 0)
  dwg.add(dwg.image(os.path.join(args.input_folder, "smr_predixcan_p4_proportion.png"), _c, _1_size))
  _label(dwg, _c, "b", (0, 50))

  dwg.save()
  t = os.path.join(args.output_folder, "supp-fig-smr-predixcan-coloc.png")
  to_png(_p, t)
  os.remove(_p)


def fig_8_deprecated_2(args):
  _p = "_fig8.svg"
  _1_size = (700, 700)
  _size = (3*_1_size[0] + 110, 2*_1_size[1] + 50)
  
  dwg = svgwrite.Drawing(_p, size=_size)
  dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))
  
  #first row
  _c = (10, 0)
  dwg.add(dwg.image(os.path.join(args.input_folder, "smr_predixcan_p3_proportion.png"), _c, _1_size))
  _label(dwg, _c, "a", (0, 50))

  _c = _advance_cursor(_c, _1_size[0]+50, 0)
  dwg.add(dwg.image(os.path.join(args.input_folder, "smr_predixcan_p3_proportion_b.png"), _c, _1_size))
  _label(dwg, _c, "b", (0, 50))

  _c = _advance_cursor(_c, _1_size[0]+50, 0)
  dwg.add(dwg.image(os.path.join(args.input_folder, "smr_h_predixcan_p3_proportion_b.png"), _c, _1_size))
  _label(dwg, _c, "c", (0, 50))

  _c = (10, _1_size[1]+50)
  dwg.add(dwg.image(os.path.join(args.input_folder, "smr_predixcan_p4_proportion.png"), _c, _1_size))
  _label(dwg, _c, "d", (0, 50))

  _c = _advance_cursor(_c, _1_size[0]+50, 0)
  dwg.add(dwg.image(os.path.join(args.input_folder, "smr_predixcan_p4_proportion_b.png"), _c, _1_size))
  _label(dwg, _c, "e", (0, 50))

  _c = _advance_cursor(_c, _1_size[0]+50, 0)
  dwg.add(dwg.image(os.path.join(args.input_folder, "smr_h_predixcan_p4_proportion_b.png"), _c, _1_size))
  _label(dwg, _c, "f", (0, 50))

  dwg.save()
  t = os.path.join(args.output_folder, "supp-fig-smr-predixcan-coloc.png")
  to_png(_p, t)
  os.remove(_p)

def fig_metaxcan_application_framework(args):
    _p = "_fig_maf.svg"
    _1_size = (2289, 1506)
    _2_size = (2193, 1098)
    _size = (_1_size[0]+100, _1_size[1] +50+ _2_size[1] )

    dwg = svgwrite.Drawing(_p, size=_size)
    dwg.add(dwg.rect(insert=(0, 0), size=_size, fill="rgb(255,255,255)"))

    _c = (100, 0)
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "metaxcan_framework.png"), _c, _1_size))
    _label(dwg, _c, "a", (-75, 100), style="font-size:100;font-family:Arial;font-weight:bold;stroke:black;stroke-width:1;fill:black")

    _c = (100+(_1_size[0]-_2_size[0])/2, _1_size[1] + 50)
    dwg.add(dwg.image(os.path.join(args.input_folder_2, "metaxcan_application.png"), _c, _2_size))
    _label(dwg, (100, _1_size[1] + 50), "b", (-75,100), style="font-size:100;font-family:Arial;font-weight:bold;stroke:black;stroke-width:1;fill:black")

    dwg.save()
    t = os.path.join(args.output_folder, "fig-metaxcan-framework-application.png")
    to_png(_p, t)
    os.remove(_p)

########################################################################################################################
def run(args):
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    fig_1(args)
    fig_2(args)
    fig_3(args)
    fig_4(args)
    fig_5(args)
    fig_6(args)
    fig_supp_coloc(args)
    fig_7_vertical(args)
    fig_metaxcan_application_framework(args)
    #fig_8(args)

if __name__ == "__main__":
    class Dummy(object):
        def __init__(self):
            self.output_folder = "results/paper_plots"
            self.input_folder = "results/plots"
            self.input_folder_2 = "data/images"

    args = Dummy()
    run(args)
