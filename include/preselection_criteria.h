#ifndef PRESELECTION_CRITERIA_H
#define PRESELECTION_CRITERIA_H

//--- trigger ---//
//diphoton trigger

//--- EE-EB gap ---//
#define CRITERION_EE_EB_GAP_LOWER             1.4442
#define CRITERION_EE_EB_GAP_UPPER             1.566

//--- photon ---//
#define CRITERION_LEADING_PHOTON_PT           35.
#define CRITERION_SUBLEADING_PHOTON_PT        25.
#define CRITERION_PHOTON_ETA                  2.5
#define CRITERION_PHOTON_IDMVA               -0.7

//--- electron ---//
// electron ID: cut based medium ID
#define CRITERION_ELECTRON_PT                 10.
#define CRITERION_ELECTRON_ETA                2.4
#define CRITERION_ELECTRON_ISO                0.2 //deltaR(e, photon)

//--- muon ---//
// muon ID: cut based tight ID
#define CRITERION_MUON_PT                     5.
#define CRITERION_MUON_ETA                    2.4
#define CRITERION_MUON_ISO                    0.2 //deltaR(mu, photon)
#define CRITERION_MUON_RELATIVE_PFISO_R04     0.25

//--- jet ---//
// jet ID: flashgg tight ID
#define CRITERION_JET_PT                      25.
#define CRITERION_JET_ETA                     2.4
#define CRITERION_JET_ISO                     0.4 //deltaR(jet, lepton/photon)

#endif
