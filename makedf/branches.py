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
    "rec.hdr.nbnbinfo",
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
    shwbranch + "end.x",   shwbranch + "end.y",   shwbranch + "end.z",
    shwbranch + "dir.x",   shwbranch + "dir.y",   shwbranch + "dir.z",
    shwbranch + 'conversion_gap', 
    shwbranch + "density",
    shwbranch + "open_angle",
    shwbranch + 'bestplane_for_energy', shwbranch + 'bestplane_for_dedx',
    shwbranch + 'bestplane_dEdx', shwbranch + 'bestplane_energy',
    shwbranch + 'plane.0.dEdx',   shwbranch + 'plane.1.dEdx', shwbranch + 'plane.2.dEdx',
    shwbranch + 'plane.0.energy', shwbranch + 'plane.1.energy', shwbranch + 'plane.2.energy',
    shwbranch + 'plane.0.nHits',  shwbranch + 'plane.1.nHits',  shwbranch + 'plane.2.nHits',
    shwbranch + "len",
    shwbranch + "truth.eff",
    shwbranch + "truth.pur",
]


presel_branches = [
    "rec.slc.is_clear_cosmic",
    "rec.slc.vertex.x", "rec.slc.vertex.y", "rec.slc.vertex.z",
    "rec.slc.nu_score",
    "rec.slc.barycenterFM.flashTime",
    "rec.slc.barycenterFM.score",
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

# making this branches separately from slcbranches since ICARUS CAFs do not have rec.slc.barycenterFM.score
barycenterFMbranches = [
    "rec.slc.barycenterFM.chargeTotal",
    "rec.slc.barycenterFM.flashTime",
    "rec.slc.barycenterFM.flashPEs",
    "rec.slc.barycenterFM.chi2",
    "rec.slc.barycenterFM.score",
    "rec.slc.barycenterFM.deltaZ",
]

correctedflashbranches = [
    "rec.slc.correctedOpFlash.OpFlashT0",
    "rec.slc.correctedOpFlash.NuToFLight",
    "rec.slc.correctedOpFlash.NuToFCharge",
    "rec.slc.correctedOpFlash.OpFlashT0Corrected",
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

mevprtltruth = "rec.mc.prtl."
mevprtltruthbranches = [
   mevprtltruth + "gen",
   mevprtltruth + "time",
   mevprtltruth + "E",
   mevprtltruth + "M",
   mevprtltruth + "flux_weight",
   mevprtltruth + "ray_weight",
   mevprtltruth + "decay_weight",
   mevprtltruth + "decay_length",
   mevprtltruth + "allowed_decay_fraction",
   mevprtltruth + "cryostat",
   mevprtltruth + "position.x",
   mevprtltruth + "position.y",
   mevprtltruth + "position.z",
   mevprtltruth + "momentum.x",
   mevprtltruth + "momentum.y",
   mevprtltruth + "momentum.z",
   mevprtltruth + "enter.x",
   mevprtltruth + "enter.y",
   mevprtltruth + "enter.z",
   mevprtltruth + "exit.x",
   mevprtltruth + "exit.y",
   mevprtltruth + "exit.z",
   mevprtltruth + "start.x",
   mevprtltruth + "start.y",
   mevprtltruth + "start.z",
   mevprtltruth + "C1",
   mevprtltruth + "C2",
   mevprtltruth + "C3",
   mevprtltruth + "C4",
   mevprtltruth + "C5",
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

# ===============  SPINE branches =============== #
# - SPINE reco interaction branches
spineint = "rec.dlp."
spineint_branches = [
    spineint + b for b in [
        "cathode_offset",                       # Distance from the cathode [cm].
        "depositions_sum",                      # Sum of deposited (de-ghosted) energy [MeV].
        "flash_hypo_pe",                        # Total PE of the hypothesized flash from OpT0Finder.
        "flash_total_pe",                       # Total PE of the matched flash (uses OpT0Finder, same order as flash_ids).
        "id",                                   # Interaction ID (dense enumeration starting from 0 within the event).
        "is_cathode_crosser",                   # Whether the interaction is a cathode-crosser (some particle crosses the cathode).
        "is_contained",                         # Whether the interaction is contained within some margin from the detector walls (see SPINE configuration)
        "is_fiducial",                          # Whether the interaction is in the fiducial volume (see SPINE configuration).
        "is_flash_matched",                     # Whether the flash is matched to the interaction (uses OpT0Finder).
        "is_matched",                           # Whether the interaction is matched to a true interaction.
        "is_time_contained",                    # Whether the particle is time-contained (within the "in-time" region of the drift window).
        "is_truth",                             # Whether the interaction is a truth interaction.
        "num_particles",                        # The number of particles in the interaction.
        "num_primary_particles",                # The number of primary particles in the interaction.
        "particle_counts.0",                    # The number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 0.
        "particle_counts.1",                    # The number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 1.
        "particle_counts.2",                    # The number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 2.
        "particle_counts.3",                    # The number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 3.
        "particle_counts.4",                    # The number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 4.
        "particle_counts.5",                    # The number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 5.
        "primary_particle_counts.0",            # The number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 0.
        "primary_particle_counts.1",            # The number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 1.
        "primary_particle_counts.2",            # The number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 2.
        "primary_particle_counts.3",            # The number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 3.
        "primary_particle_counts.4",            # The number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 4.
        "primary_particle_counts.5",            # The number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 5.
        "size",                                 # The size of the interaction (number of voxels).
#        "topology",                             # Topology of the interaction (e.g. "0g0e1mu0pi2p") considering only primaries.
        "vertex.0",                             # Vertex of the interaction in detector coordinates (X coordinate).
        "vertex.1",                             # Vertex of the interaction in detector coordinates (Y coordinate).
        "vertex.2",                             # Vertex of the interaction in detector coordinates (Z coordinate).
#        "particles",                            # Particles in the interaction.
    ]
]

spineint_flashids_branches = [
    spineint + "flash_ids"                      # Flash IDs for the matched flashes (uses OpT0Finder, lowest score first).
]

spineint_flashscores_branches = [
    spineint + "flash_scores"                   # Score of the matched flashes (uses OpT0Finder, same order as flash_ids).  
]

spineint_flashtimes_branches = [
    spineint + "flash_times"                    # Time of the matched flashes (uses OpT0Finder, same order as flash_ids).
]

spineint_flashvolumeids_branches = [
    spineint + "flash_volume_ids"               # Volume IDs of the matched flashes (uses OpT0Finder, same order as flash_ids).
]

spineint_matched_branches = [
    spineint + "match_ids"                      # Interaction IDs of the considered matches (correspond to true interaction IDs).
]

spineint_matchovrlp_branches = [
    spineint + "match_overlaps"                 # Intersection over union (IoU) of space points of the considered matches.
]

spineint_moduleids_branches = [
    spineint + "module_ids"                     # Module IDs (cryostat) of the interaction.
]

spineint_particleids_branches = [
    spineint + "particle_ids"                   # Particle IDs in the interaction (see 'id' attribute of SRParticleDLP).
]

spineint_primaryparticleids_branches = [
    spineint + "primary_particle_ids"           # Primary particle IDs in the interaction (see 'id' attribute of SRParticleDLP).
]

# - SPINE true interaction branches
spinetint = "rec.dlp_true."

spinetint_branches = [
    spinetint + b for b in [
        "bjorken_x",                            # Bjorken x of the neutrino interaction.
        "cathode_offset",                       # Distance from the cathode [cm].
#        "creation_process",                     # Creation process associated with the neutrino (see MCParticle).
        "current_type",                         # Current type of the neutrino (see MCNeutrino; 0=CC, 1=NC).
        "depositions_adapt_q_sum",              # Total tagged (reco non-ghost) charge deposited [ADC].
        "depositions_adapt_sum",                # Total tagged (reco non-ghost) energy deposited [MeV].
        "depositions_g4_sum",                   # Total energy deposited energy at the G4 level [MeV].
        "depositions_q_sum",                    # Total tagged (true non-ghost) charge deposited [ADC].
        "depositions_sum",                      # Total tagged (true non-ghost) energy deposited [MeV].
        "distance_travel",                      # Distance traveled by the neutrino from production to the interaction.
        "energy_init",                          # Initial energy of the neutrino.
        "energy_transfer",                      # Energy transfer (Q0) of the neutrino interaction.
        "flash_hypo_pe",                        # Total PE of the hypothesized flash.
        "flash_total_pe",                       # Total PE of the matched flash.
        "hadronic_invariant_mass",              # Hadronic invariant mass of the neutrino.
        "id",                                   # Interaction ID (dense enumeration starting from 0 within the event).
        "inelasticity",                         # Inelasticity of the neutrino interaction.
        "interaction_mode",                     # Interaction mode of the neutrino (see MCNeutrino).
        "interaction_type",                     # Interaction type of the neutrino (see MCNeutrino).
        "is_cathode_crosser",                   # Whether the interaction is a cathode-crosser (some particle crosses the cathode).
        "is_contained",                         # Whether the interaction is contained within some margin from the detector walls.
        "is_fiducial",                          # Whether the interaction is in the fiducial volume (see SPINE configuration).
        "is_flash_matched",                     # Whether the flash is matched to the interaction.
        "is_matched",                           # Whether the interaction is matched to a true interaction.
        "is_time_contained",                    # Whether the particle is time-contained (within the "in-time" region of the drift window).
        "is_truth",                             # Whether the interaction is a truth interaction.
        "lepton_p",                             # Momentum of the lepton in the interaction.
        "lepton_pdg_code",                      # PDG code of the lepton in the interaction.
        "lepton_track_id",                      # Track ID of the lepton in the neutrino interaction.
        "mct_index",                            # Index of the neutrino in the original MCTruth array.
        "momentum.0",                           # Momentum (vector) of the neutrino (X coordinate).
        "momentum.1",                           # Momentum (vector) of the neutrino (Y coordinate).
        "momentum.2",                           # Momentum (vector) of the neutrino (Z coordinate).
        "momentum_transfer",                    # Momentum transfer (Q^2) of the neutrino interaction.
        "momentum_transfer_mag",                # Momentum transfer (Q3) of the neutrino interaction.
        "nu_id",                                # Neutrino ID (-1 = not a neutrino, 0 = first neutrino, 1 = second neutrino, etc.).
        "nucleon",                              # Hit nucleon in the neutrino interaction (see MCNeutrino).
        "num_particles",                        # Number of particles in the interaction.
        "num_primary_particles",                # Number of primary particles in the interaction.
        "orig_id",                              # Original ID of the interaction (from Geant4).
        "particle_counts.0",                    # Number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 0.
        "particle_counts.1",                    # Number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 1.
        "particle_counts.2",                    # Number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 2.
        "particle_counts.3",                    # Number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 3.
        "particle_counts.4",                    # Number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 4.
        "particle_counts.5",                    # Number of particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 5.
        "pdg_code",                             # PDG code of the neutrino.
        "position.0",                           # Position of the neutrino interaction (X coordinate).
        "position.1",                           # Position of the neutrino interaction (Y coordinate).
        "position.2",                           # Position of the neutrino interaction (Z coordinate).
        "primary_particle_counts.0",            # Number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 0.
        "primary_particle_counts.1",            # Number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 1.
        "primary_particle_counts.2",            # Number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 2.
        "primary_particle_counts.3",            # Number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 3.
        "primary_particle_counts.4",            # Number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 4.
        "primary_particle_counts.5",            # Number of primary particles of each type in the interaction (photon, electron, muon, pion, proton, kaon), type 5.
        "quark",                                # Hit quark in the neutrino interaction (see MCNeutrino).
        "reco_vertex.0",                        # Vertex of the interaction in detector coordinates (reco) (X coordinate).
        "reco_vertex.1",                        # Vertex of the interaction in detector coordinates (reco) (Y coordinate).
        "reco_vertex.2",                        # Vertex of the interaction in detector coordinates (reco) (Z coordinate).
        "size",                                 # Number of true non-ghost true-tagged space points.
        "size_adapt",                           # Number of reco non-ghost true-tagged space points.
        "size_g4",                              # Number of (rasterized) g4 energy depositions (no detector effects).
        "target",                               # Target in the neutrino interaction.
        "theta",                                # Angle between the neutrino and the outgoing lepton [radians].
#        "topology",                             # Topology of the interaction (e.g. "0g0e1mu0pi2p") considering only primaries.
        "track_id",                             # Track ID of the neutrino interaction.
        "vertex.0",                             # Vertex of the interaction in detector coordinates (truth) (X coordinate).
        "vertex.1",                             # Vertex of the interaction in detector coordinates (truth) (Y coordinate).
        "vertex.2",                             # Vertex of the interaction in detector coordinates (truth) (Z coordinate).
#        "particles",                            # Particles in the interaction.
    ]
]

spinetint_flashids_branches = [
    spinetint + "flash_ids"                      # Flash IDs for the matched flashes (uses OpT0Finder, lowest score first).
]
spinetint_flashscores_branches = [
    spinetint + "flash_scores"                   # Flash score for the matched flashes (uses OpT0Finder, same order as flash_ids).  
]
spinetint_flashtimes_branches = [
    spinetint + "flash_times"                    # Time of the matched flashes (uses OpT0Finder, same order as flash_ids).
]
spinetint_flashvolumeids_branches = [
    spinetint + "flash_volume_ids"               # Volume IDs of the matched flashes (uses OpT0Finder, same order as flash_ids).
]
spinetint_matched_branches = [
    spinetint + "match_ids"                      # Interaction IDs of the considered matches (correspond to reco interaction IDs).
]
spinetint_matchovrlp_branches = [
    spinetint + "match_overlaps"                 # Intersection over union (IoU) of space points of the considered matches.
]
spinetint_moduleids_branches = [
    spinetint + "module_ids"                     # Module IDs (cryostat) of the interaction.
]
spinetint_particleids_branches = [
    spinetint + "particle_ids"                   # Particle IDs in the interaction (see 'id' attribute of SRParticleTruthDLP).
]
spinetint_primaryparticleids_branches = [
    spinetint + "primary_particle_ids"           # Primary particle IDs in the interaction (see 'id' attribute of SRParticleTruthDLP).
]


# - SPINE reco particle branches
spinepart = "rec.dlp.particles."

spinepart_branches = [
    spinepart + b for b in [
        "axial_spread",                         # Axial spread of the particle.
        "calo_ke",                              # Calorimetric kinetic energy.
        "cathode_offset",                       # Distance from the cathode [cm].
        "chi2_per_pid.0",                       # Chi2 score for PID hypothesis 0.
        "chi2_per_pid.1",                       # Chi2 score for PID hypothesis 1.
        "chi2_per_pid.2",                       # Chi2 score for PID hypothesis 2.
        "chi2_per_pid.3",                       # Chi2 score for PID hypothesis 3.
        "chi2_per_pid.4",                       # Chi2 score for PID hypothesis 4.
        "chi2_per_pid.5",                       # Chi2 score for PID hypothesis 5.
        "chi2_pid",                             # PID from the chi2-based PID.
        "csda_ke",                              # Continuous-slowing-down-approximation kinetic energy.
        "csda_ke_per_pid.0",                    # CSDA kinetic energy for PID 0.
        "csda_ke_per_pid.1",                    # CSDA kinetic energy for PID 1.
        "csda_ke_per_pid.2",                    # CSDA kinetic energy for PID 2.
        "csda_ke_per_pid.3",                    # CSDA kinetic energy for PID 3.
        "csda_ke_per_pid.4",                    # CSDA kinetic energy for PID 4.
        "csda_ke_per_pid.5",                    # CSDA kinetic energy for PID 5.
        "depositions_sum",                      # Sum of deposited (de-ghosted) energy [MeV].
        "directional_spread",                   # Directional spread of the particle.
        "end_dir.0",                            # Unit direction vector calculated at the particle end point (X coordinate).
        "end_dir.1",                            # Unit direction vector calculated at the particle end point (Y coordinate).
        "end_dir.2",                            # Unit direction vector calculated at the particle end point (Z coordinate).
        "end_point.0",                          # End point (vector) of the particle (X coordinate).
        "end_point.1",                          # End point (vector) of the particle (Y coordinate).
        "end_point.2",                          # End point (vector) of the particle (Z coordinate).
        "id",                                   # Particle ID (dense enumeration starting from 0 within the event).
        "interaction_id",                       # Parent interaction ID.
        "is_cathode_crosser",                   # Whether the particle is a cathode-crosser.
        "is_contained",                         # Whether the particle is contained within some margin from the detector walls (see SPINE configuration).
        "is_matched",                           # Whether the particle is matched.
        "is_primary",                           # Whether the particle is a primary particle.
        "is_time_contained",                    # Whether the particle is time-contained (within the "in-time" region of the drift window).
        "is_truth",                             # Whether the particle is a truth particle.
        "is_valid",                             # Whether the particle passes thresholds and counts towards the topology.
        "ke",                                   # Kinetic energy from best energy estimator (CSDA, calorimetric, or MCS).
        "length",                               # Length of the particle.
        "mass",                                 # Mass of the particle (according to assigned PID).
        "mcs_ke",                               # Multiple Coulomb scattering kinetic energy.
        "mcs_ke_per_pid.0",                     # MCS kinetic energy for PID 0.
        "mcs_ke_per_pid.1",                     # MCS kinetic energy for PID 1.
        "mcs_ke_per_pid.2",                     # MCS kinetic energy for PID 2.
        "mcs_ke_per_pid.3",                     # MCS kinetic energy for PID 3.
        "mcs_ke_per_pid.4",                     # MCS kinetic energy for PID 4.
        "mcs_ke_per_pid.5",                     # MCS kinetic energy for PID 5.
        "momentum.0",                           # Momentum (vector) of the particle (X coordinate).
        "momentum.1",                           # Momentum (vector) of the particle (Y coordinate).
        "momentum.2",                           # Momentum (vector) of the particle (Z coordinate).
        "num_fragments",                        # Number of fragments comprising the particle.
        "p",                                    # Momentum magnitude.
        "pdg_code",                             # PDG code of the particle.
        "pid",                                  # Particle ID (see Pid_t enumeration).
        "pid_scores.0",                         # PID 0 softmax scores.
        "pid_scores.1",                         # PID 1 softmax scores.
        "pid_scores.2",                         # PID 2 softmax scores.
        "pid_scores.3",                         # PID 3 softmax scores.
        "pid_scores.4",                         # PID 4 softmax scores.
        "pid_scores.5",                         # PID 5 softmax scores.
        "primary_scores.0",                     # Primary score for PID 0.
        "primary_scores.1",                     # Primary score for PID 1.
        "shape",                                # Semantic type of the particle (see Shape_t enumeration).
        "size",                                 # The size of the particle (number of voxels).
        "start_dedx",                           # dE/dx at the start of the particle.
        "start_dir.0",                          # Unit direction vector calculated at the particle start point (X coordinate).
        "start_dir.1",                          # Unit direction vector calculated at the particle start point (Y coordinate).
        "start_dir.2",                          # Unit direction vector calculated at the particle start point (Z coordinate).
        "start_point.0",                        # Start point (vector) of the particle (X coordinate).
        "start_point.1",                        # Start point (vector) of the particle (Y coordinate).
        "start_point.2",                        # Start point (vector) of the particle (Z coordinate).
        "start_straightness",                   # Straightness at the start of the particle.
        "vertex_distance"                       # Distance from the vertex.
    ]
]

spinepart_fragmentids_branches = [
    spinepart + "fragment_ids"                  # Fragment IDs comprising the particle.
]

spinepart_matched_branches = [
    spinepart + "match_ids"                     # Particle IDs of the considered matches (correspond to true particle IDs).
]

spinepart_matchovrlp_branches = [
    spinepart + "match_overlaps"                # Intersection over union (IoU) of space points of the considered matches.
]

spinepart_moduleids_branches = [
    spinepart + "module_ids"                    # Module IDs (cryostat) of the particle.
]

spinepart_ppnids_branches = [
    spinepart + "ppn_ids"                       # PPN IDs of the particle.
]

# - SPINE true particle branches
spinetpart = "rec.dlp_true.particles."

spinetpart_branches = [
    spinetpart + b for b in [
#        "ancestor_creation_process",           # Geant4 creation process of the ancestor particle.
        "ancestor_pdg_code",                    # PDG code of the ancestor particle.
        "ancestor_position.0",                  # Position of the ancestor particle (X coordinate).
        "ancestor_position.1",                  # Position of the ancestor particle (Y coordinate).
        "ancestor_position.2",                  # Position of the ancestor particle (Z coordinate).
        "ancestor_t",                           # Time of the ancestor particle.
        "ancestor_track_id",                    # Track ID of the ancestor particle.
        "calo_ke",                              # Calorimetric kinetic energy.
        "cathode_offset",                       # Distance from the cathode [cm].
#        "children_counts",                     # Number of children of the particle.
#        "creation_process",                    # Geant4 creation process associated with the particle (see MCParticle).
        "csda_ke",                              # Continuous-slowing-down-approximation kinetic energy.
        "csda_ke_per_pid.0",                    # CSDA kinetic energy per PID (PID 0).
        "csda_ke_per_pid.1",                    # CSDA kinetic energy per PID (PID 1).
        "csda_ke_per_pid.2",                    # CSDA kinetic energy per PID (PID 2).
        "csda_ke_per_pid.3",                    # CSDA kinetic energy per PID (PID 3).
        "csda_ke_per_pid.4",                    # CSDA kinetic energy per PID (PID 4).
        "csda_ke_per_pid.5",                    # CSDA kinetic energy per PID (PID 5).
        "depositions_adapt_q_sum",              # Total tagged (reco non-ghost) charge deposited [ADC].
        "depositions_adapt_sum",                # Total tagged (reco non-ghost) energy deposited [MeV].
        "depositions_g4_sum",                   # Total energy deposited energy at the G4 level [MeV].
        "depositions_q_sum",                    # Total tagged (true non-ghost) charge deposited [ADC].
        "depositions_sum",                      # Total tagged (true non-ghost) energy deposited [MeV].
        "distance_travel",                      # Distance traveled by the neutrino from production to the interaction.
        "end_dir.0",                            # Unit direction vector at the particle end point (X coordinate).
        "end_dir.1",                            # Unit direction vector at the particle end point (Y coordinate).
        "end_dir.2",                            # Unit direction vector at the particle end point (Z coordinate).
        "end_momentum.0",                       # Momentum vector at the particle end (X coordinate).
        "end_momentum.1",                       # Momentum vector at the particle end (Y coordinate).
        "end_momentum.2",                       # Momentum vector at the particle end (Z coordinate).
        "end_p",                                # Momentum magnitude of the particle at the end.
        "end_point.0",                          # End point vector of the particle (X coordinate).
        "end_point.1",                          # End point vector of the particle (Y coordinate).
        "end_point.2",                          # End point vector of the particle (Z coordinate).
        "end_position.0",                       # End position vector of the particle (X coordinate).
        "end_position.1",                       # End position vector of the particle (Y coordinate).
        "end_position.2",                       # End position vector of the particle (Z coordinate).
        "end_t",                                # Time of the particle at the end.
        "energy_deposit",                       # Energy deposited by the particle.
        "energy_init",                          # Initial energy of the particle.
        "first_step.0",                         # X coordinate of the first step of the particle.
        "first_step.1",                         # Y coordinate of the first step oftheparticle.
        "first_step.2",                         # Z coordinate of the first stepoftheparticle.
        "group_id",                             # Group ID of the particle.
        "group_primary",                        # Whether the particle is a primary within its group.
        "id",                                   # Particle ID (dense enumeration starting from 0 within the event).
        "interaction_id",                       # Parent interaction ID.
        "interaction_primary",                  # Whether the particle is a primary within its interaction (equivalent to is_primary).
        "is_cathode_crosser",                   # Whether the particle is a cathode-crosser.
        "is_contained",                         # Whether the particle is contained within some margin from the detector walls (see SPINE configuration).
        "is_matched",                           # Whether the particle is matched.
        "is_primary",                           # Whether the particle is a primary particle.
        "is_time_contained",                    # Whether the particle is time-contained (within the "in-time" region of the drift window).
        "is_truth",                             # Whether the particle is a truth particle.
        "is_valid",                             # Whether the particle passes thresholds and counts towards topology.
        "ke",                                   # Kinetic energy from best energy estimator (CSDA, calorimetric, or MCS).
        "last_step.0",                          # X coordinate of the last step of the particle.
        "last_step.1",                          # Y coordinate of the coordinates of the last step of the particle.
        "last_step.2",                          # Z coordinate of the coordinates of the last step of the particle.
        "length",                               # Length of the particle.
        "mass",                                 # Mass of the particle.
        "mcs_ke",                               # Multiple Coulomb scattering kinetic energy.
        "mcs_ke_per_pid.0",                     # MCS kinetic energy per PID.
        "mcs_ke_per_pid.1",                     # MCS kinetic energy per PID.
        "mcs_ke_per_pid.2",                     # MCS kinetic energy per PID.
        "mcs_ke_per_pid.3",                     # MCS kinetic energy per PID.
        "mcs_ke_per_pid.4",                     # MCS kinetic energy per PID.
        "mcs_ke_per_pid.5",                     # MCS kinetic energy per PID.
        "mcst_index",                           # MCST index.
        "mct_index",                            # Index of the particle in the original MCTruth array.
        "momentum.0",                           # Mmomentum (vector) of the particle (X coordinate).
        "momentum.1",                           # Momentum (vector) of the particle (Y coordinate).
        "momentum.2",                           # Momentum (vector) of the particle (Z coordinate).
        "nu_id",                                # Neutrino ID (-1 = not a neutrino, 0 = first neutrino, 1 = second neutrino, etc.).
        "num_fragments",                        # Number of particle fragments.
        "num_voxels",                           # Number of voxels comprising the particle.
        "orig_group_id",                        # Original group ID of the particle.
        "orig_id",                              # Original ID of the particle.
        "orig_interaction_id",                  # Interaction ID as it was stored in the parent LArCV file under the interaction_id attribute.
        "orig_parent_id",                       # Parent ID as it was stored in the parent LArCV file under the parent_id attribute.
        "p",                                    # Momentum magnitude.
#        "parent_creation_process",             # Geant4 creation process of the parent particle.
        "parent_id",                            # Parent particle ID.
        "parent_pdg_code",                      # PDG code of the parent particle.
        "parent_position.0",                    # Position of the parent particle (X coordinate).
        "parent_position.1",                    # Position of the parent particle (Y coordinate).
        "parent_position.2",                    # Position of the parent particle (Z coordinate).
        "parent_t",                             # Time of the parent particle.
        "parent_track_id",                      # Track ID of the parent particle.
        "pdg_code",                             # PDG code of the particle.
        "pid",                                  # Particle ID (see Pid enumeration).
        "position.0",                           # Position of the particle (X coordinate).
        "position.1",                           # Position of the particle (Y coordinate).
        "position.2",                           # Position of the particle (Z coordinate).
        "reco_end_dir.0",                       # Direction vector at the reconstructed end point of the particle (X coordinate).
        "reco_end_dir.1",                       # Direction vector at the reconstructed end point of the particle (Y coordinate).
        "reco_end_dir.2",                       # Direction vector at the reconstructed end point of the particle (Z coordinate).
        "reco_ke",                              # Kinetic energy estimator using reconstructed information.
        "reco_length",                          # Reconstructed length of the particle.
        "reco_momentum.0",                      # Reconstructed momentum (vector) of the particle.
        "reco_momentum.1",                      # Reconstructed momentum (vector) of the particle.
        "reco_momentum.2",                      # Reconstructed momentum (vector) of the particle.
        "reco_start_dir.0",                     # Direction vector at the reconstructed start point of the particle (X coordinate).
        "reco_start_dir.1",                     # Direction vector at the reconstructed start point of the particle (Y coordinate).
        "reco_start_dir.2",                     # Direction vector at the reconstructed start point of the particle (Z coordinate).
        "shape",                                # Semantic type of the particle (see SemanticType enumeration).
        "size",                                 # Size of the particle (number of voxels).
        "size_adapt",                           # Size of the particle using adaptive thresholds (number of voxels).
        "size_g4",                              # Size of the particle at the G4 level (number of voxels).
        "start_dir.0",                          # Unit direction vector at the particle start point (X coordinate).
        "start_dir.1",                          # Unit direction vector at the particle start point (Y coordinate).
        "start_dir.2",                          # Unit direction vector at the particle start point (Z coordinate).
        "start_point.0",                        # Position of the start point (vector) of the particle (X coordinate).
        "start_point.1",                        # Position of the start point (vector) of the particle (Y coordinate).
        "start_point.2",                        # Position of the start point (vector) of the particle (Z coordinate).
        "t",                                    # Time of the particle.
        "track_id"                              # Track ID of the particle.
    ]   
]   


spinetpart_childrenid_branches = [
    spinetpart + "children_id"                  # List of particle ID of children particles.
]

spinetpart_fragmentids_branches = [
    spinetpart + "fragment_ids"                 # Fragment IDs comprising the particle.
]

spinetpart_matched_branches = [
    spinetpart + "match_ids"                    # Particle IDs of the considered matches (correspond to reco particle IDs).
]

spinetpart_matchovrlp_branches = [
    spinetpart + "match_overlaps"               # Intersection over union (IoU) of space points of the considered matches.
]

spinetpart_moduleids_branches = [
    spinetpart + "module_ids"                   # Module IDs (cryostat) of the particle.
]

spinetpart_originchildrenid_branches = [
    spinetpart + "orig_children_id"             # Original ID of the children particles.
]
