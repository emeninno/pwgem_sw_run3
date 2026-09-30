import datetime
import sys
sys.path.append("../common/");
import numpy as np
import pandas as pd
import ROOT
from ROOT import TFile, TCanvas, TPad, TLegend, TPaveText
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kOrange, kMagenta, kCyan
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kOpenCircle, kOpenSquare, kOpenTriangleUp, kOpenDiamond, kOpenCross, kFullStar, kOpenStar, kOpenTriangleDown, kFullDiamond, kFullCross
from painter import make_common_style
from histo_manager import get_ratio
from signal_extractor import get_significance
gStyle.SetPalette(55);
gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);

#__________________________________________________________
def draw_mee_uls_ls(filename, taskname, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    print(filename)
    print(taskname)
    print(type(rootdire))
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);
    h1m_uls.Scale(1, "width");
    h1m_bkg.Scale(1, "width");
    h1m_sig.Scale(1, "width");

    h1m_uls.SetDirectory(0);
    h1m_bkg.SetDirectory(0);
    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_uls, False);
    ROOT.SetOwnership(h1m_bkg, False);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_uls, 20, 1.2, kBlack, 2, 0);
    make_common_style(h1m_bkg, 25, 1.2, kBlue+1, 2, 0);
    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = h1m_uls.GetMaximum() * 50;
    y_range_min = max(h1m_sig.GetMinimum() * 0.2, 2e-8);

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{ee}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1m_uls.Draw("E0same");
    h1m_bkg.Draw("E0same");
    h1m_sig.Draw("E0same");

    leg = TLegend(0.65, 0.8, 0.85, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    leg.AddEntry(h1m_uls, "ULS = S + B", "P");
    leg.AddEntry(h1m_bkg, "B = 2R #sqrt{#it{N}_{++} #it{N}_{#minus#minus}}", "P");
    leg.AddEntry(h1m_sig, "S = ULS #minus B", "P");
    leg.Draw("");
    ROOT.SetOwnership(leg,False);

    txt = TPaveText(0.2, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE");
    txt.AddText("PbPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs("{0}_raw_dndm_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_uls_ls_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
def draw_mee_sbratio(filename, taskname, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);


    h1sb = get_ratio(h1m_sig, h1m_bkg, "G");
    make_common_style(h1sb, 20, 1.2, kRed+1, 2, 0);
    h1sb.SetDirectory(0);
    ROOT.SetOwnership(h1sb, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    #y_range_max = h1sb.GetMaximum() * 10;
    #y_range_min = h1sb.GetMinimum() * 0.1;
    y_range_max = 2e+3;
    y_range_min = 1e-2;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("S/B");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1sb.Draw("E0same");

    txt = TPaveText(0.2, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE");
    txt.AddText("PbPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_mee_sbratio_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    #c1.SaveAs("{0}_mee_sbratio_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sbratio_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_significance(filename, taskname, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    rootdire.ls();
    hs_uls = rootdire.Get("hs_uls");
    hs_bkg = rootdire.Get("hs_bkg");
    hs_sig = rootdire.Get("hs_sig");
    h1z = None;
    for key in rootdire.GetListOfKeys():
    	if key.GetName() == "Event/after/hZvtx":
        	h1z = key.ReadObj()
        	break

    if not h1z:
    	print("Error: Could not find object 'Event/after/hZvtx'.")
    	return
    #h1z = rootdire.FindObject("Event/after/hZvtx");
    #if not h1z:
    #	print("Error: Could not find object 'Event/after/hZvtx'.")
    #	return
    nev = h1z.GetEntries();

    bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
    bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
    hs_uls.GetAxis(1).SetRange(bin0, bin1);
    hs_bkg.GetAxis(1).SetRange(bin0, bin1);
    hs_sig.GetAxis(1).SetRange(bin0, bin1);

    bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_uls.GetAxis(2).SetRange(bin0, bin1);
    hs_bkg.GetAxis(2).SetRange(bin0, bin1);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);

    h1m_uls = hs_uls.Projection(0);
    h1m_bkg = hs_bkg.Projection(0);
    h1m_sig = hs_sig.Projection(0);
    h1m_bkg.Scale(nev);
    h1m_sig.Scale(nev);

    h1significance = get_significance(h1m_sig, h1m_bkg);
    make_common_style(h1significance, 20, 1.2, kRed+1, 2, 0);
    h1significance.SetDirectory(0);
    ROOT.SetOwnership(h1significance, False);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = 1e+4;
    y_range_min = 1e+0;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("significance = S/#sqrt{S + 2B}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    h1significance.Draw("E0same");

    txt = TPaveText(0.2, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE");
    txt.AddText("PbPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_mee_significance_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    #c1.SaveAs("{0}_mee_significance_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_significance_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_sbratio_multiple(filename, tasknames, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = 2e+3;
    y_range_min = 1e-2;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("S/B");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.65, 0.7, 0.85, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    ROOT.SetOwnership(leg,False);

    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
        bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
        hs_uls.GetAxis(2).SetRange(bin0, bin1);
        hs_bkg.GetAxis(2).SetRange(bin0, bin1);
        hs_sig.GetAxis(2).SetRange(bin0, bin1);

        h1m_uls = hs_uls.Projection(0);
        h1m_bkg = hs_bkg.Projection(0);
        h1m_sig = hs_sig.Projection(0);

        h1sb = get_ratio(h1m_sig, h1m_bkg, "G");
        h1sb.SetDirectory(0);
        ROOT.SetOwnership(h1sb, False);
        make_common_style(h1sb, 20, 1.2, colors[it], 2, 0);

        h1sb.Draw("E0same");

        taskname_tmp = taskname.replace("dielectron-qc_", "")
        if taskname == "dielectron-qc":
            taskname_tmp = "default";

        leg.AddEntry(h1sb, taskname_tmp, "P");

    txt = TPaveText(0.2, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE");
    txt.AddText("PbPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_mee_sbratio_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    #c1.SaveAs("{0}_mee_sbratio_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_sbratio_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_significance_multiple(filename, tasknames, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = 1e+4;
    y_range_min = 1e+0;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("significance = S/#sqrt{S + 2B}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.65, 0.7, 0.85, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    ROOT.SetOwnership(leg,False);

    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");
        h1z = rootdire.FindObject("Event/after/hZvtx");
        nev = h1z.GetEntries();

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
        bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
        hs_uls.GetAxis(2).SetRange(bin0, bin1);
        hs_bkg.GetAxis(2).SetRange(bin0, bin1);
        hs_sig.GetAxis(2).SetRange(bin0, bin1);

        h1m_uls = hs_uls.Projection(0);
        h1m_bkg = hs_bkg.Projection(0);
        h1m_sig = hs_sig.Projection(0);
        h1m_bkg.Scale(nev);
        h1m_sig.Scale(nev);

        h1significance = get_significance(h1m_sig, h1m_bkg);
        h1significance.SetDirectory(0);
        ROOT.SetOwnership(h1significance, False);
        make_common_style(h1significance, 20, 1.2, colors[it], 2, 0);

        h1significance.Draw("E0same");

        taskname_tmp = taskname.replace("dielectron-qc_", "")
        if taskname == "dielectron-qc":
            taskname_tmp = "default";
        leg.AddEntry(h1significance, taskname_tmp, "L");

    txt = TPaveText(0.2, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE");
    txt.AddText("PbPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_mee_significance_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    #c1.SaveAs("{0}_mee_significance_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_mee_significance_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_mee_multiple(filename, tasknames, ptmin, ptmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    colors = [kBlack, kRed+1, kGreen+2, kBlue+1, kYellow+1, kOrange+1, kMagenta+1, kCyan+1];

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.16, 0.02, 0.1, 0.02);
    c1.SetTicks(1,1);
    c1.SetLogy(1);

    y_range_max = 0;
    y_range_min = 0;

    frame1 = c1.DrawFrame(0., y_range_min, 4, y_range_max);
    frame1.GetXaxis().SetTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{m}_{ee}} (GeV/#it{c}^{2})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.8);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);

    leg = TLegend(0.65, 0.7, 0.85, 0.95);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.035);
    ROOT.SetOwnership(leg,False);

    ntask = len(tasknames);
    for it in range(0, ntask):
        taskname = tasknames[it];
        rootdire = rootfile.Get(taskname);
        hs_uls = rootdire.Get("hs_uls");
        hs_bkg = rootdire.Get("hs_bkg");
        hs_sig = rootdire.Get("hs_sig");

        bin0 = hs_uls.GetAxis(1).FindBin(ptmin + 1e-3);
        bin1 = hs_uls.GetAxis(1).FindBin(ptmax - 1e-3);
        hs_uls.GetAxis(1).SetRange(bin0, bin1);
        hs_bkg.GetAxis(1).SetRange(bin0, bin1);
        hs_sig.GetAxis(1).SetRange(bin0, bin1);

        bin0 = hs_uls.GetAxis(2).FindBin(dcamin + 1e-3);
        bin1 = hs_uls.GetAxis(2).FindBin(dcamax - 1e-3);
        hs_uls.GetAxis(2).SetRange(bin0, bin1);
        hs_bkg.GetAxis(2).SetRange(bin0, bin1);
        hs_sig.GetAxis(2).SetRange(bin0, bin1);

        h1m_uls = hs_uls.Projection(0);
        h1m_bkg = hs_bkg.Projection(0);
        h1m_sig = hs_sig.Projection(0);
        h1m_uls.Scale(1, "width");
        h1m_bkg.Scale(1, "width");
        h1m_sig.Scale(1, "width");

        if it == 0:
            y_range_max = h1m_uls.GetMaximum() * 50;
            y_range_min = h1m_sig.GetMinimum() * 0.2;

        make_common_style(h1m_sig, 20, 1.2, colors[it], 2, 0);
        h1m_sig.Draw("E0same");

        taskname_tmp = taskname.replace("dielectron-qc_", "")
        if taskname == "dielectron-qc":
            taskname_tmp = "default";

        leg.AddEntry(h1m_sig, taskname_tmp, "P");

    frame1.GetYaxis().SetRangeUser(y_range_min, y_range_max);
    txt = TPaveText(0.2, 0.65, 0.5, 0.95, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE");
    txt.AddText("PbPb at #sqrt{#it{s}_{NN}} = 5.36 TeV");
    txt.AddText("#it{p}_{T,e} > 0.2 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/#it{{c}}".format(ptmin, ptmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    leg.Draw("");

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_raw_dndm_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    #c1.SaveAs("{0}_raw_dndm_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, ptmin, ptmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndm_multiple_pt{1:2.1f}_{2:2.1f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, ptmin, ptmax, dcamin, dcamax, suffix));

    rootfile.Close();

#__________________________________________________________
def draw_ptee_upc(filename, taskname, cen1, cen2, mmin, mmax, dcamin, dcamax, suffix=""):
    rootfile = TFile(filename, "READ");
    rootdire = rootfile.Get(taskname);
    hs_sig = rootdire.Get("hs_sig");

    bin0 = hs_sig.GetAxis(0).FindBin(mmin + 1e-3);
    bin1 = hs_sig.GetAxis(0).FindBin(mmax - 1e-3);
    hs_sig.GetAxis(0).SetRange(bin0, bin1);

    bin0 = hs_sig.GetAxis(2).FindBin(dcamin + 1e-3);
    bin1 = hs_sig.GetAxis(2).FindBin(dcamax - 1e-3);
    hs_sig.GetAxis(2).SetRange(bin0, bin1);
    h1m_sig = hs_sig.Projection(1);
    h1m_sig.Scale(1, "width");

    h1m_sig.SetDirectory(0);
    ROOT.SetOwnership(h1m_sig, False);

    make_common_style(h1m_sig, 20, 1.2, kRed+1, 2, 0);

    c1 = TCanvas("c0", "c0", 0, 0, 800, 800);
    c1.SetMargin(0.18, 0.02, 0.1, 0.05);
    c1.SetTicks(1,1);
    #c1.SetLogy(1);

    y_range_max = h1m_sig.GetMaximum() * 5;
    y_range_min = 0.0;

    frame1 = c1.DrawFrame(0., y_range_min, 1, y_range_max);
    frame1.GetXaxis().SetTitle("#it{p}_{T,ee} (GeV/#it{c})");
    frame1.GetYaxis().SetTitle("#frac{1}{#it{N}_{ev}} #frac{d#it{N}}{d#it{p}_{T,ee}} (GeV/#it{c})^{#minus1}");
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(2.1);
    frame1.GetXaxis().SetTitleSize(0.04);
    frame1.GetYaxis().SetTitleSize(0.04);
    frame1.GetXaxis().SetLabelSize(0.04);
    frame1.GetYaxis().SetLabelSize(0.04);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    frame1.GetXaxis().SetMoreLogLabels(1);
    ROOT.SetOwnership(frame1,False);
    h1m_sig.Draw("E0same");

    txt = TPaveText(0.22, 0.62, 0.5, 0.92, "NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(12);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.035);
    txt.AddText("ALICE");
    txt.AddText("{0}#minus{1}% Pb#minusPb at #sqrt{{#it{{s}}_{{NN}}}} = 5.36 TeV".format(cen1, cen2));
    txt.AddText("#it{p}_{T,e} > 0.4 GeV/#it{c}, |#it{#eta}_{e}| < 0.8");
    txt.AddText("|#it{y}_{ee}| < 0.8");
    txt.AddText("{0:3.2f} < #it{{m}}_{{ee}} < {1:3.2f} GeV/#it{{c}}^{{2}}".format(mmin, mmax));
    txt.AddText("{0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma".format(dcamin, dcamax));
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    date = datetime.date.today().strftime("%Y%m%d");
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    #c1.SaveAs("{0}_raw_dndpt_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.eps".format(date, mmin, mmax, dcamin, dcamax, suffix));
    #c1.SaveAs("{0}_raw_dndpt_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.pdf".format(date, mmin, mmax, dcamin, dcamax, suffix));
    c1.SaveAs("{0}_raw_dndpt_m{1:3.2f}_{2:3.2f}GeV_dca3d{3:2.1f}_{4:2.1f}sigma{5}.png".format(date, mmin, mmax, dcamin, dcamax, suffix));

    rootfile.Close();
#__________________________________________________________
#__________________________________________________________
#__________________________________________________________
if __name__ == "__main__":
    filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC25_pass1_0090.root";
    taskname = "dielectron_pre_occi_0_5000_cent_0090_detadphiCut";
    ptmin = 0;
    ptmax = 10;
    dcamin = 0;
    dcamax = 10;
    suffix = "";
    #suffix = "_0_5000_cent_0090_detadphiCut";
    draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0,  0.5, suffix);
    draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 2.0, 10.0, suffix);

    draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    draw_mee_sbratio(filename, taskname, 0.0, 10.0, 0.0,  0.5, suffix);
    draw_mee_sbratio(filename, taskname, 0.0, 10.0, 2.0, 10.0, suffix);

    draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    draw_mee_significance(filename, taskname, 0.0, 10.0, 0.0,  0.5, suffix);
    draw_mee_significance(filename, taskname, 0.0, 10.0, 2.0, 10.0, suffix);

    #tasknames = [
    #    "dielectron-qc"
    #    ,"dielectron-qc_itsibany"
    #    ,"dielectron-qc_newsel8"
    #];
    #draw_mee_sbratio_multiple(filename, tasknames, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_significance_multiple(filename, tasknames, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_multiple(filename, tasknames, 0.0, 10.0, 0.0, 10.0, suffix);


    ## for peripheral PbPb
    #filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_5070.root";
    #taskname = "dielectron-qc_5070_TOFreq_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_5070_TOFreq_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 1.1, 2.5, 0.0,  0.5, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 2.7, 3.2, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 50, 70, 2.7, 3.2, 0.0,  0.5, suffix);

    #filename = "mee_ptee_dcaee_PbPb_5.36TeV_LHC23zzh_pass3_7090.root";
    #taskname = "dielectron-qc_7090_TOFreq_newsel8";
    #suffix = "_PbPb_5.36TeV_LHC23zzh_pass3_7090_TOFreq_newsel8";
    #draw_mee_uls_ls(filename, taskname, 0.0, 10.0, 0.0, 10.0, suffix);
    #draw_mee_uls_ls(filename, taskname, 0.0, 0.1, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 1.1, 2.5, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 1.1, 2.5, 0.0,  0.5, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 2.7, 3.2, 0.0, 10.0, suffix);
    #draw_ptee_upc(filename, taskname, 70, 90, 2.7, 3.2, 0.0,  0.5, suffix);


#__________________________________________________________
