import os
import sys
import shutil
import numpy as np
import yaml
import argparse
import ROOT
import subprocess

SEGFAULT_PROCESS_RETURNCODE = -11


segfaulting_code = "import ctypes ; ctypes.string_at(0)"  # https://codegolf.stackexchange.com/a/4694/115779
try:
    subprocess.run(["python3", "-c", segfaulting_code],
                   check=True)
except subprocess.CalledProcessError as err:
    if err.returncode == SEGFAULT_PROCESS_RETURNCODE:
        print("probably segfaulted")
    else:
        print(f"crashed for other reasons: {err.returncode}")
else:
    print("ok")

ROOT.gROOT.SetBatch(True)
from ROOT import TFile, TDirectory, TH1F, TH2F, THnSparseF
from dilepton_analyzer import DileptonAnalyzer

parser = argparse.ArgumentParser('Example program');
parser.add_argument("-i", "--input" , default="AnalysisResults.root", type=str, help="path to the root file you want to analyze", required=True)
parser.add_argument("-c", "--config", default="config.yml", type=str, help="path to the *.yml configuration file", required=True)
parser.add_argument("-t", "--type"  , default="data" , type=str, help="run type [data or mc]", required=True)
parser.add_argument("-s", "--suffix"  , default="" , type=str, help="suffix for output file name", required=False)
args = parser.parse_args();

config = "";
with open(args.config, "r", encoding="utf-8") as config_yml:
    config = yaml.safe_load(config_yml)
#print(config);
system = config["common"]["system"];
energy = config["common"]["energy"]; #float
suffix = args.suffix;

period = "";
is_mc = False;
if "mc" in args.type:
    is_mc = True;
    period = config["period_mc"];
else:
    is_mc = False;
    period = config["period_data"];
pass_number = config["pass_number"];
#__________________________________________________________________
def set_range_00(hs):
    hs.GetAxis(0).SetRange(0, 0);
    hs.GetAxis(1).SetRange(0, 0);
    hs.GetAxis(2).SetRange(0, 0);
#__________________________________________________________________
def run_mee_ptee_dcaee(filename, tasknames, arr_mee, arr_ptee, arr_dcaee):
    print(sys._getframe().f_code.co_name);
    print("filename = ", filename);
    print("arr_mee = ", arr_mee);
    print("arr_ptee = ", arr_ptee);
    print("arr_dcaee = ", arr_dcaee);
    __delta = 1e-3;

    rootfile = TFile.Open(filename, "READ");
    outname = "";
    if "pp" in system:
        outname = "mee_ptee_dcaee_{0}_{1:3.1f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
    else:
        outname = "mee_ptee_dcaee_{0}_{1:3.2f}TeV_{2}_{3}{4}.root".format(system, energy, period, pass_number, suffix)
    print("outname = ", outname);
    outfile = TFile(outname, "RECREATE");

    #for itask in range(0, 1):
    for itask in range(0, len(tasknames)):
        taskname = tasknames[itask];
        print("analyzing ... ", taskname); 

        rootdir_task = rootfile.Get(taskname);
        if rootdir_task == None:
            print(taskname, "does not exist. skip.");
            continue;

        rootdir_event = rootdir_task.Get("Event");
        rootdir_pair = rootdir_task.Get("Pair");

        rootdir_event = rootdir_task.Get("Event/after");
        rootdir_pair = rootdir_task.Get("Pair");
        #rootdir_event.ls();
        #rootdir_pair .ls();
        outdir = outfile.mkdir(taskname);

        h1z = rootdir_event.Get("hZvtx");
        h1centrality = rootdir_event.Get("hCentFT0C");
        outdir.WriteTObject(h1z);
        outdir.WriteTObject(h1centrality);
        nev = h1z.GetEntries();

        hs_uls_same  = rootdir_pair.Get("same/uls/hs") .Clone("hs_uls_same")
        hs_lspp_same = rootdir_pair.Get("same/lspp/hs").Clone("hs_lspp_same")
        hs_lsmm_same = rootdir_pair.Get("same/lsmm/hs").Clone("hs_lsmm_same")
        hs_uls_mix   = rootdir_pair.Get("mix/uls/hs")  .Clone("hs_uls_mix")
        hs_lspp_mix  = rootdir_pair.Get("mix/lspp/hs") .Clone("hs_lspp_mix")
        hs_lsmm_mix  = rootdir_pair.Get("mix/lsmm/hs") .Clone("hs_lsmm_mix")
        hs_uls_same .Sumw2();
        hs_lspp_same.Sumw2();
        hs_lsmm_same.Sumw2();
        hs_uls_mix  .Sumw2();
        hs_lspp_mix .Sumw2();
        hs_lsmm_mix .Sumw2();

        #outdir.WriteTObject(hs_uls_same);
        #outdir.WriteTObject(hs_lspp_same);
        #outdir.WriteTObject(hs_lsmm_same);
        #outdir.WriteTObject(hs_uls_mix);
        #outdir.WriteTObject(hs_lspp_mix);
        #outdir.WriteTObject(hs_lsmm_mix);

        set_range_00(hs_uls_same );
        set_range_00(hs_lspp_same);
        set_range_00(hs_lsmm_same);
        set_range_00(hs_uls_mix  );
        set_range_00(hs_lspp_mix );
        set_range_00(hs_lsmm_mix );

        __ndim = 3;
        nbin = np.array([len(arr_mee)-1, len(arr_ptee)-1, len(arr_dcaee)-1], dtype=np.int32);
        xmin = np.array([arr_mee[0], arr_ptee[0], arr_dcaee[0]], dtype=np.float64);
        xmax = np.array([arr_mee[-1], arr_ptee[-1], arr_dcaee[-1]], dtype=np.float64);
        
        # Debug prints for nbin, xmin, xmax
        #print(f"nbin: {nbin}")
        #print(f"xmin: {xmin}")
        #print(f"xmax: {xmax}")
        hs_uls = THnSparseF("hs_uls", "#it{N}_{ee}^{all}/#it{N}_{ev};#it{m}_{ee} (GeV/#it{c}^{2});#it{p}_{T,ee} (GeV/#it{c});DCA_{ee}^{3D} (#sigma);", __ndim, nbin, xmin, xmax);
        # Check histogram initialization
        if not hs_uls:
            raise RuntimeError("Failed to initialize hs_uls histogram.")
        hs_uls.Sumw2();

        hs_uls.SetBinEdges(0, arr_mee);
        hs_uls.SetBinEdges(1, arr_ptee);
        hs_uls.SetBinEdges(2, arr_dcaee);

        # Debug prints for bin edges
        #print(f"arr_mee: {arr_mee}")
        #print(f"arr_ptee: {arr_ptee}")
        #print(f"arr_dcaee: {arr_dcaee}")
        #print(f"hs_uls.GetAxis(0).GetXmin(): {hs_uls.GetAxis(0).GetXmin()}, hs_uls.GetAxis(0).GetXmax(): {hs_uls.GetAxis(0).GetXmax()}")
        #print(f"hs_uls.GetAxis(1).GetXmin(): {hs_uls.GetAxis(1).GetXmin()}, hs_uls.GetAxis(1).GetXmax(): {hs_uls.GetAxis(1).GetXmax()}")
        #print(f"hs_uls.GetAxis(2).GetXmin(): {hs_uls.GetAxis(2).GetXmin()}, hs_uls.GetAxis(2).GetXmax(): {hs_uls.GetAxis(2).GetXmax()}")
        hs_R = hs_uls.Clone("hs_R");
        hs_R.SetTitle("R factor");
        hs_bkg = hs_uls.Clone("hs_bkg");
        hs_bkg.SetTitle("#it{N}_{ee}^{bkg}/#it{N}_{ev}");
        hs_sig = hs_uls.Clone("hs_sig");
        hs_sig.SetTitle("#it{N}_{ee}^{signal}/#it{N}_{ev}");

        for idca in range(0, len(arr_dcaee)-1):
            dca_min = arr_dcaee[idca];
            dca_max = arr_dcaee[idca+1];
            dca_center = (dca_min + dca_max)/2;
            bin_dca1 = hs_uls_same.GetAxis(2).FindBin(dca_min + __delta);
            bin_dca2 = hs_uls_same.GetAxis(2).FindBin(dca_max - __delta);
            #print("bin_dca1: ", bin_dca1)
            #print("bin_dca2: ", bin_dca2)

            hs_uls_same .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lspp_same.GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lsmm_same.GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_uls_mix  .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lspp_mix .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            hs_lsmm_mix .GetAxis(2).SetRange(bin_dca1, bin_dca2);
            #h2R = TH2F("h2R_dca{0}".format(idca), "R factor in {0:2.1f} < DCA_{{ee}}^{{3D}} < {1:2.1f} #sigma;#it{{m}}_{{ee}} (GeV/#it{{c}}^{{2}});#it{{p}}_{{T,ee}} (GeV/#it{{c}})".format(dca_min, dca_max) , len(arr_mee)-1, arr_mee, len(arr_ptee)-1, arr_ptee);
            #h2R.Sumw2();
            #h2R.SetContour(100);

            for ipt in range(0, len(arr_ptee)-1):
                pt_min = arr_ptee[ipt];
                pt_max = arr_ptee[ipt+1];
                pt_center = (pt_min + pt_max)/2;
                bin_pt1 = hs_uls_same.GetAxis(1).FindBin(pt_min + __delta);
                bin_pt2 = hs_uls_same.GetAxis(1).FindBin(pt_max - __delta);
                #print("bin_pt1: ", bin_pt1)
                #print("bin_pt2: ", bin_pt2)
                hs_uls_same .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lspp_same.GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lsmm_same.GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_uls_mix  .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lspp_mix .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                hs_lsmm_mix .GetAxis(1).SetRange(bin_pt1, bin_pt2);
                h1m_uls_same  = hs_uls_same.Projection(0);
                h1m_lspp_same = hs_lspp_same.Projection(0);
                h1m_lsmm_same = hs_lsmm_same.Projection(0);
                h1m_uls_mix   = hs_uls_mix.Projection(0);
                h1m_lspp_mix  = hs_lspp_mix.Projection(0);
                h1m_lsmm_mix  = hs_lsmm_mix.Projection(0);
                h1m_uls_same .SetName("h1m_uls_same_pt{0}_dca{1}".format(ipt, idca));
                h1m_lspp_same.SetName("h1m_lspp_same_pt{0}_dca{1}".format(ipt, idca));
                h1m_lsmm_same.SetName("h1m_lsmm_same_pt{0}_dca{1}".format(ipt, idca));
                h1m_uls_mix  .SetName("h1m_uls_mix_pt{0}_dca{1}".format(ipt, idca));
                h1m_lspp_mix .SetName("h1m_lspp_mix_pt{0}_dca{1}".format(ipt, idca));
                h1m_lsmm_mix .SetName("h1m_lsmm_mix_pt{0}_dca{1}".format(ipt, idca));

                h1m_uls_same .SetTitle("#it{{N}}_{{+-}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_same.SetTitle("#it{{N}}_{{++}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_same.SetTitle("#it{{N}}_{{--}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_uls_mix  .SetTitle("#it{{N}}_{{+-}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_mix .SetTitle("#it{{N}}_{{++}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_mix .SetTitle("#it{{N}}_{{--}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));

                h1m_uls_same .SetXTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
                h1m_lspp_same.SetXTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
                h1m_lsmm_same.SetXTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
                h1m_uls_mix  .SetXTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
                h1m_lspp_mix .SetXTitle("#it{m}_{ee} (GeV/#it{c}^{2})");
                h1m_lsmm_mix .SetXTitle("#it{m}_{ee} (GeV/#it{c}^{2})");

                h1m_uls_same .SetYTitle("counts per bin");
                h1m_lspp_same.SetYTitle("counts per bin");
                h1m_lsmm_same.SetYTitle("counts per bin");
                h1m_uls_mix  .SetYTitle("counts per bin");
                h1m_lspp_mix .SetYTitle("counts per bin");
                h1m_lsmm_mix .SetYTitle("counts per bin");

                ana = DileptonAnalyzer(h1m_uls_same, h1m_lspp_same, h1m_lsmm_same, h1m_uls_mix, h1m_lspp_mix, h1m_lsmm_mix, arr_mee);
                [h1m_uls_same_rebin, h1m_lspp_same_rebin, h1m_lsmm_same_rebin, h1m_uls_mix_rebin, h1m_lspp_mix_rebin, h1m_lsmm_mix_rebin, h1R, h1bkg, h1sig] = ana.run();    
                h1m_uls_same_rebin .SetName("h1m_uls_same_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lspp_same_rebin.SetName("h1m_lspp_same_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lsmm_same_rebin.SetName("h1m_lsmm_same_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_uls_mix_rebin  .SetName("h1m_uls_mix_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lspp_mix_rebin .SetName("h1m_lspp_mix_pt{0}_dca{1}_rebin".format(ipt, idca));
                h1m_lsmm_mix_rebin .SetName("h1m_lsmm_mix_pt{0}_dca{1}_rebin".format(ipt, idca));

                h1m_uls_same_rebin .SetTitle("#it{{N}}_{{+-}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_same_rebin.SetTitle("#it{{N}}_{{++}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_same_rebin.SetTitle("#it{{N}}_{{--}}^{{same}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_uls_mix_rebin  .SetTitle("#it{{N}}_{{+-}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lspp_mix_rebin .SetTitle("#it{{N}}_{{++}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1m_lsmm_mix_rebin .SetTitle("#it{{N}}_{{--}}^{{mix}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));

                h1R  .SetName("h1R_pt{0}_dca{1}"  .format(ipt, idca));
                h1bkg.SetName("h1bkg_pt{0}_dca{1}".format(ipt, idca));
                h1sig.SetName("h1sig_pt{0}_dca{1}".format(ipt, idca));
                h1R  .SetTitle("R factor in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1bkg.SetTitle("#it{{N}}_{{ee}}^{{bkg}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));
                h1sig.SetTitle("#it{{N}}_{{ee}}^{{signal}} in {0:2.1f} < #it{{p}}_{{T,ee}} < {1:2.1f} GeV/c, {2:2.1f} < DCA_{{ee}}^{{3D}} < {3:2.1f} #sigma".format(pt_min, pt_max, dca_min, dca_max));

                h1R  .SetYTitle("#it{R} = #frac{#it{N}_{+#minus}^{mix}}{2 #sqrt{#it{N}_{++}^{mix} #it{N}_{#minus#minus}^{mix}}}");
                h1bkg.SetYTitle("counts per bin");
                h1sig.SetYTitle("counts per bin");

                #outdir.WriteTObject(h1m_uls_same_rebin );
                #outdir.WriteTObject(h1m_lspp_same_rebin);
                #outdir.WriteTObject(h1m_lsmm_same_rebin);
                #outdir.WriteTObject(h1m_uls_mix_rebin  );
                #outdir.WriteTObject(h1m_lspp_mix_rebin );
                #outdir.WriteTObject(h1m_lsmm_mix_rebin );
                #outdir.WriteTObject(h1R);
                #outdir.WriteTObject(h1bkg);
                #outdir.WriteTObject(h1s);
                for im in range(0, len(arr_mee)-1):
                    m_min = arr_mee[im];
                    m_max = arr_mee[im+1];
                    m_center = (m_min + m_max)/2;
                    global_bin_id = hs_uls.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_R  .GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_bkg.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    global_bin_id = hs_sig.GetBin(np.array([m_center, pt_center, dca_center], dtype=np.float64), True);
                    #print(f"global_bin_id: {global_bin_id}")
                    #print(f"m_center: {m_center}, pt_center: {pt_center}, dca_center: {dca_center}")
                    if global_bin_id == 0:
                        print(f"Invalid global_bin_id for coordinates: m_center={m_center}, pt_center={pt_center}, dca_center={dca_center}")
                        continue  # Skip this iteration
                    hs_uls.SetBinContent(global_bin_id, h1m_uls_same_rebin.GetBinContent(im+1));
                    hs_uls.SetBinError(global_bin_id, h1m_uls_same_rebin.GetBinError(im+1));
                    hs_R  .SetBinContent(global_bin_id, h1R.GetBinContent(im+1));
                    hs_R  .SetBinError(global_bin_id, h1R.GetBinError(im+1));
                    hs_bkg.SetBinContent(global_bin_id, h1bkg.GetBinContent(im+1));
                    hs_bkg.SetBinError(global_bin_id, h1bkg.GetBinError(im+1));
                    hs_sig.SetBinContent(global_bin_id, h1sig.GetBinContent(im+1));
                    hs_sig.SetBinError(global_bin_id, h1sig.GetBinError(im+1));
                    #h2R.SetBinContent(im+1, ipt+1, h1R.GetBinContent(im+1));
                    #h2R.SetBinError(im+1, ipt+1, h1R.GetBinError(im+1));
                    #print("End of setting bin content and bin error for hs histograms")
            #outdir.WriteTObject(h2R);
            #print("End of the dca loop");
        hs_uls.Scale(1/nev);
        hs_bkg.Scale(1/nev);
        hs_sig.Scale(1/nev);
        outdir.WriteTObject(hs_uls);
        outdir.WriteTObject(hs_R  );
        outdir.WriteTObject(hs_bkg);
        outdir.WriteTObject(hs_sig);

    outfile.Close();
    rootfile.Close();
#__________________________________________________________________
if __name__ == "__main__":
    filename = args.input;

    arr_mee = np.array(config["common"]["mee_bin"],dtype=float);
    arr_ptee = np.array(config["common"]["ptee_bin"],dtype=float);
    arr_dcaee = np.array(config["common"]["dcaee_bin"],dtype=float);

    if is_mc :
        tasknames = config["mc"]["tasknames"]
        #print(tasknames);
    else:
        tasknames = config["data"]["tasknames"]
        #print(tasknames);
        run_mee_ptee_dcaee(filename, tasknames, arr_mee, arr_ptee, arr_dcaee);


