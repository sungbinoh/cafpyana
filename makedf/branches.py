mchdrbranches = [
    "rec.hdr.pot",
    "rec.hdr.first_in_subrun",
    "rec.hdr.ismc",
    "rec.hdr.run",
    "rec.hdr.subrun",
    "rec.hdr.ngenevt",
    "rec.hdr.evt",
    "rec.hdr.proc",
    "rec.hdr.cluster",
    "rec.hdr.fno",
]

hdrbranches = [
    "rec.hdr.pot",
    "rec.hdr.first_in_subrun",
    "rec.hdr.ismc",
    "rec.hdr.run",
    "rec.hdr.subrun",
    "rec.hdr.ngenevt",
    "rec.hdr.evt",
    "rec.hdr.proc",
    "rec.hdr.cluster",
    "rec.hdr.fno",
    "rec.hdr.noffbeambnb",
]

trigger_info_branches = [
    "rec.hdr.triggerinfo.trigger_id",
    "rec.hdr.triggerinfo.gate_id",
    "rec.hdr.triggerinfo.trigger_count",
    "rec.hdr.triggerinfo.gate_count",
    "rec.hdr.triggerinfo.gate_delta",
    "rec.hdr.triggerinfo.global_trigger_time",
    "rec.hdr.triggerinfo.prev_global_trigger_time",
    "rec.hdr.triggerinfo.global_trigger_det_time",
    "rec.hdr.triggerinfo.trigger_within_gate",
    "rec.hdr.triggerinfo.beam_gate_det_time",
    "rec.hdr.triggerinfo.beam_gate_time_abs",
]

opflashbranches = [
    "rec.opflashes.firsttime",
    "rec.opflashes.time",
    "rec.opflashes.totalpe",
]

numipotbranches = [
    "rec.hdr.numiinfo.spill_time_s",
    "rec.hdr.numiinfo.spill_time_ns",
    "rec.hdr.numiinfo.TRTGTD",
    "rec.hdr.numiinfo.TORTGT",
    "rec.hdr.numiinfo.daq_gates",
]

bnbpotbranches = [
    "rec.hdr.bnbinfo.TOR860",
    "rec.hdr.bnbinfo.TOR875",
]

sbndframebranches = [
    "rec.sbnd_frames.frameApplyAtCaf",
    "rec.sbnd_frames.frameHltBeamGate",
    "rec.sbnd_frames.frameHltCrtt1",
    "rec.sbnd_frames.frameTdcBes",
    "rec.sbnd_frames.frameTdcCrtt1",
    "rec.sbnd_frames.frameTdcRwm",
    "rec.sbnd_frames.timingType",
]

sbndtimingbranches = [
    "rec.sbnd_timings.hltBeamGate",
    "rec.sbnd_timings.hltCrtt1",
    "rec.sbnd_timings.hltEtrig",
    "rec.sbnd_timings.rawDAQHeaderTimestamp",
    "rec.sbnd_timings.tdcBes",
    "rec.sbnd_timings.tdcCrtt1",
    "rec.sbnd_timings.tdcEtrig",
    "rec.sbnd_timings.tdcRwm",
]

trueparticlenames = [
    "start_process",
    "end_process",
    "pdg",
    "startE",
    "endE",
    "start.x", "start.y", "start.z",
    "end.x", "end.y", "end.z",
    "genp.x", "genp.y", "genp.z",
    "length",
    "G4ID",
    "parent",
    "cont_tpc",
    "genE",
    "interaction_id",
    "crosses_tpc",
]

trueparticlebranches = [
    *["rec.true_particles.%s" % s for s in trueparticlenames],
    "rec.true_particles.plane.0.2.visE",
    "rec.true_particles.plane.1.2.visE",
]

crtvetobranches = [      
    "rec.sbnd_crt_veto.V0",
    "rec.sbnd_crt_veto.V1",
    "rec.sbnd_crt_veto.V2",
    "rec.sbnd_crt_veto.V3",
    "rec.sbnd_crt_veto.V4"            
]

crtspbranches = [      
    "rec.crt_spacepoints.pe",
    "rec.crt_spacepoints.position.x",
    "rec.crt_spacepoints.position.y",
    "rec.crt_spacepoints.position.z",
    "rec.crt_spacepoints.position_err.x",
    "rec.crt_spacepoints.position_err.y",
    "rec.crt_spacepoints.position_err.z",
    "rec.crt_spacepoints.time",
    "rec.crt_spacepoints.time_err"
]

crthitbranches = [
  "rec.crt_hits.time",
  "rec.crt_hits.t1",
  "rec.crt_hits.t0",
  "rec.crt_hits.pe",
  "rec.crt_hits.plane",
]


pfpbranch = "rec.slc.reco.pfp."
trkbranch = pfpbranch + "trk."
shwbranch = pfpbranch + "shw."

pfobranches = [
    pfpbranch + "pfochar.chgendfrac",
    pfpbranch + "pfochar.chgfracspread",
    pfpbranch + "pfochar.linfitdiff",
    pfpbranch + "pfochar.linfitlen",
    pfpbranch + "pfochar.linfitgaplen",
    pfpbranch + "pfochar.linfitrms",
    pfpbranch + "pfochar.openanglediff",
    pfpbranch + "pfochar.pca2ratio",
    pfpbranch + "pfochar.pca3ratio", 
    pfpbranch + "pfochar.vtxdist" 
]

pfpbranches = [
    pfpbranch + "parent_is_primary",
    pfpbranch + "slcID",
    pfpbranch + "trackScore",
    pfpbranch + "parent",
    pfpbranch + "id",
    pfpbranch + "t0",
] + pfobranches

pfp_daughter_branch = [
    pfpbranch + "daughters"
]

trkbranches = [
    trkbranch + "producer",
    trkbranch + "start.x", trkbranch + "start.y", trkbranch + "start.z",
    trkbranch + "end.x", trkbranch + "end.y", trkbranch + "end.z",
    trkbranch + "dir.x", trkbranch + "dir.y", trkbranch + "dir.z",
    trkbranch + "len",
    trkbranch + "rangeP.p_muon",
    trkbranch + "mcsP.fwdP_muon",
    trkbranch + "rangeP.p_pion",
    trkbranch + "mcsP.fwdP_pion",
    trkbranch + "rangeP.p_proton",
    trkbranch + "mcsP.fwdP_proton",
    trkbranch + "bestplane",
    trkbranch + "crthit.distance",
    trkbranch + "crthit.hit.time",
    trkbranch + "crthit.hit.pe",
    trkbranch + "chi2pid.2.pid_ndof",
    trkbranch + "chi2pid.2.chi2_muon",
    trkbranch + "chi2pid.2.chi2_proton",
    trkbranch + "chi2pid.2.pida",
] + pfpbranches

trkmcsbranches = [
  trkbranch + "mcsP.seg_length",
  trkbranch + "mcsP.seg_scatter_angles",
]

shwbranches = [
    shwbranch + "start.x", shwbranch + "start.y", shwbranch + "start.z",
    shwbranch + "end.x",   shwbranch + "end.y", shwbranch + "end.z",
    shwbranch + 'conversion_gap', 
    shwbranch + "density",
    shwbranch + "open_angle",
    shwbranch + 'bestplane',
    shwbranch + 'bestplane_dEdx', shwbranch + 'bestplane_energy',
    shwbranch + 'plane.0.dEdx',   shwbranch + 'plane.1.dEdx', shwbranch + 'plane.2.dEdx',
    shwbranch + 'plane.0.energy', shwbranch + 'plane.1.energy', shwbranch + 'plane.2.energy',
    shwbranch + 'plane.0.nHits',  shwbranch + 'plane.1.nHits',  shwbranch + 'plane.2.nHits',
    shwbranch + "len",
    shwbranch + "truth.eff",
    shwbranch + "truth.pur",
]

trkhitadcbranches = [
  trkbranch + "calo.2.points.adcs"
]

trkhitbranches_perplane = lambda IPLANE : [
    trkbranch + "calo.%i.points.dedx"% IPLANE,
    trkbranch + "calo.%i.points.dqdx"% IPLANE,
    trkbranch + "calo.%i.points.pitch"% IPLANE,
    trkbranch + "calo.%i.points.integral"% IPLANE,
    trkbranch + "calo.%i.points.rr"% IPLANE,
    trkbranch + "calo.%i.points.phi"% IPLANE,
    trkbranch + "calo.%i.points.efield"% IPLANE,
    trkbranch + "calo.%i.points.wire"% IPLANE,
    trkbranch + "calo.%i.points.tpc"% IPLANE,
    trkbranch + "calo.%i.points.sumadc"% IPLANE,
    trkbranch + "calo.%i.points.t"% IPLANE,
    trkbranch + "calo.%i.points.x"% IPLANE,
    trkbranch + "calo.%i.points.y"% IPLANE,
    trkbranch + "calo.%i.points.z"% IPLANE,

    #trkbranch + "calo.%i.points.width"% IPLANE,
    #trkbranch + "calo.%i.points.mult"% IPLANE,
    #trkbranch + "calo.%i.points.tdc0"% IPLANE,
]

trktruehitbranches_perplane = lambda IPLANE : [
    trkbranch + "calo.%i.points.truth.h_e"% IPLANE,
    trkbranch + "calo.%i.points.truth.h_nelec"% IPLANE,
]

trkhitbranches = trkhitbranches_perplane(2)
trkhitbranches_P1 = trkhitbranches_perplane(1)
trkhitbranches_P0 = trkhitbranches_perplane(0)

trktruehitbranches = trktruehitbranches_perplane(2)
trktruehitbranches_P1 = trktruehitbranches_perplane(1)
trktruehitbranches_P0 = trktruehitbranches_perplane(0)

#### ICARUS flat.caf does not have efield and phi for each hit so far,
trkhitbranches_perplane_icarus = lambda IPLANE : [
    trkbranch + "calo.%i.points.dedx"% IPLANE,
    trkbranch + "calo.%i.points.dqdx"% IPLANE,
    trkbranch + "calo.%i.points.pitch"% IPLANE,
    trkbranch + "calo.%i.points.integral"% IPLANE,
    trkbranch + "calo.%i.points.rr"% IPLANE,
    trkbranch + "calo.%i.points.wire"% IPLANE,
    trkbranch + "calo.%i.points.tpc"% IPLANE,
    trkbranch + "calo.%i.points.sumadc"% IPLANE,
    trkbranch + "calo.%i.points.t"% IPLANE,
    trkbranch + "calo.%i.points.x"% IPLANE,
    trkbranch + "calo.%i.points.y"% IPLANE,
    trkbranch + "calo.%i.points.z"% IPLANE,
]

trkhitbranches_icarus = trkhitbranches_perplane_icarus(2)
trkhitbranches_P1_icarus = trkhitbranches_perplane_icarus(1)
trkhitbranches_P0_icarus = trkhitbranches_perplane_icarus(0)

for n in trueparticlenames: trkbranches.append(trkbranch + "truth.p." + n)
for n in trueparticlenames: shwbranches.append(shwbranch + "truth.p." + n)

slcbranches = [
    "rec.slc.is_clear_cosmic",
    "rec.slc.vertex.x", "rec.slc.vertex.y", "rec.slc.vertex.z",
    "rec.slc.self",
    "rec.slc.tmatch.eff",
    "rec.slc.tmatch.pur",
    "rec.slc.tmatch.index",
    "rec.slc.producer",
    "rec.slc.nuid.crlongtrkdiry",
    "rec.slc.nu_score",
    "rec.slc.opt0.score",
    "rec.slc.opt0.time",
    "rec.slc.opt0.measPE",
    "rec.slc.opt0.hypoPE",
]

## making this branches separately to support sbnd flat.caf versions older than v10_10_*
barycenterFMbranches = [
    "rec.slc.barycenterFM.chargeTotal",
    "rec.slc.barycenterFM.flashTime",
    "rec.slc.barycenterFM.flashPEs",
    "rec.slc.barycenterFM.chi2",
    "rec.slc.barycenterFM.score",
]

mcbranches = [
    "rec.mc.nu.E",
    "rec.mc.nu.baseline",
    "rec.mc.nu.time",
    "rec.mc.nu.bjorkenX",
    "rec.mc.nu.inelasticityY",
    "rec.mc.nu.Q2",
    "rec.mc.nu.w",
    "rec.mc.nu.momentum.x",
    "rec.mc.nu.momentum.y",
    "rec.mc.nu.momentum.z",
    "rec.mc.nu.position.x",
    "rec.mc.nu.position.y",
    "rec.mc.nu.position.z",
    "rec.mc.nu.pdg",
    "rec.mc.nu.iscc",
    "rec.mc.nu.genie_mode",
    "rec.mc.nu.parent_pdg",
    "rec.mc.nu.parent_dcy_E",
    "rec.mc.nu.genie_evtrec_idx",
]

mcprimbranches = [
    "rec.mc.nu.prim.genE",
    "rec.mc.nu.prim.length",
    "rec.mc.nu.prim.pdg",
    "rec.mc.nu.prim.genp.x",
    "rec.mc.nu.prim.genp.y",
    "rec.mc.nu.prim.genp.z",
    "rec.mc.nu.prim.start.x", "rec.mc.nu.prim.start.y", "rec.mc.nu.prim.start.z",
    "rec.mc.nu.prim.end.x", "rec.mc.nu.prim.end.y", "rec.mc.nu.prim.end.z",
    "rec.mc.nu.prim.interaction_id", "rec.mc.nu.prim.crosses_tpc", "rec.mc.nu.prim.contained",
]

mcprimvisEbranches = [
    "rec.mc.nu.prim.plane.0.2.visE", "rec.mc.nu.prim.plane.1.2.visE"
]

mcprimdaughtersbranches = [
    "rec.mc.nu.prim.daughters",
]

slc_mcbranches = ["rec.slc.truth." + ".".join(s.split(".")[3:]) for s in mcbranches]
slc_mcprimbranches = ["rec.slc.truth." + ".".join(s.split(".")[3:]) for s in mcprimbranches]

mchbranches = [
  "rec.mc.prtl.time",
  "rec.mc.prtl.E",
  "rec.mc.prtl.M",
  "rec.mc.prtl.start.x", "rec.mc.prtl.start.y", "rec.mc.prtl.start.z",
  "rec.mc.prtl.enter.x", "rec.mc.prtl.enter.y", "rec.mc.prtl.enter.z",
  "rec.mc.prtl.exit.x", "rec.mc.prtl.exit.y", "rec.mc.prtl.exit.z",
  "rec.mc.prtl.decay_length",
  "rec.mc.prtl.allowed_decay_fraction",
  "rec.mc.prtl.C1",
  "rec.mc.prtl.C2",
  "rec.mc.prtl.C3",
  "rec.mc.prtl.C4",
  "rec.mc.prtl.C5",
]

stubbranches = [
    "rec.slc.reco.stub.vtx.x",
    "rec.slc.reco.stub.vtx.y",
    "rec.slc.reco.stub.vtx.z",
    "rec.slc.reco.stub.end.x",
    "rec.slc.reco.stub.end.y",
    "rec.slc.reco.stub.end.z",

    "rec.slc.reco.stub.efield_vtx",
    "rec.slc.reco.stub.efield_end",


    "rec.slc.reco.stub.truth.p.pdg",
    "rec.slc.reco.stub.truth.p.genE",
    "rec.slc.reco.stub.truth.p.interaction_id",
]

stubplanebranches = [
    "rec.slc.reco.stub.planes.p",
    "rec.slc.reco.stub.planes.hit_w",
    "rec.slc.reco.stub.planes.vtx_w",
    "rec.slc.reco.stub.planes.pitch",
    "rec.slc.reco.stub.planes.trkpitch",
]

stubhitbranches = [
    "rec.slc.reco.stub.planes.hits.charge",
    "rec.slc.reco.stub.planes.hits.ontrack",
    "rec.slc.reco.stub.planes.hits.wire",
]

eslc = "rec.dlp."

eslcbranches = [
    eslc + "is_neutrino",
    eslc + "nu_id",
    eslc + "num_particles",
    eslc + "num_primaries",
    eslc + "vertex.0",
    eslc + "vertex.1",
    eslc + "vertex.2",
]

eslcmatchedbranches = [
    eslc + "match",
]

eslcmatchovrlpbranches = [
    eslc + "match_overlap",
]

etruthint = "rec.dlp_true."

etruthintbranches = [
    etruthint + "id"
]

epart = "rec.dlp.particles."

eparticlebranches = [
    epart + "end_point.0",
    epart + "end_point.1",
    epart + "end_point.2",
    epart + "is_contained",
    epart + "is_primary",
    epart + "is_principal_match",
    epart + "is_valid",
    epart + "length",
    epart + "csda_ke",
    epart + "ke",
    epart + "momentum.0",
    epart + "momentum.1",
    epart + "momentum.2",
    epart + "pid",
    epart + "pid_scores.0",
    epart + "pid_scores.1",
    epart + "pid_scores.2",
    epart + "pid_scores.3",
    epart + "pid_scores.4",
    epart + "start_point.0",
    epart + "start_point.1",
    epart + "start_point.2",
    epart + "start_dir.0",
    epart + "start_dir.1",
    epart + "start_dir.2",
]

eparticlematchedbranches = [
    epart + "match",
]

eparticlematchovrlpbranches = [
    epart + "match_overlap",
]

etruthpart = "rec.dlp_true.particles."

etrueparticlebranches = [
    etruthpart + "track_id",
    etruthpart + "id",
]

etruthint = "rec.dlp_true."

etruthintbranches = [
    etruthint + "id",
    etruthint + "nu_id"
]


