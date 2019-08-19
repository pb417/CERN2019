import ROOT as r
import numpy as np
from AnalysisHelpers.ROOTUtils import DefineTree, SaveTree

#Import  the data
ch =r.TChain('DecayTree')
ch.Add("2016DstrD2KsKKLLFixed_PostTrigger0.root")

#Function GetMass out of a TLorentzVector

def GetMass (px, py, pz, E):
    p = r.TLorentzVector()
    p.SetPxPyPzE(px, py, pz, E)
    M = p.M()
    return M

def GetMDzPi(ev,daus):
    mdz = 1864.84 #MeV
    mpi = 139.57#MeV
    vecdz = r.TVector3(0,0,0)
    for v, tv in daus.iteritems():
        vecdz += r.TVector3(ev.__getattr__("Dst_DTF_D0_%s_PX"%tv),
                            ev.__getattr__("Dst_DTF_D0_%s_PY"%tv),
                            ev.__getattr__("Dst_DTF_D0_%s_PZ"%tv))
    vecpi = r.TVector3(ev.Dst_DTF_pis_PX,
                       ev.Dst_DTF_pis_PY,
                       ev.Dst_DTF_pis_PZ)
    mdz2, mpi2 = mdz*mdz, mpi*mpi
    mdzpi = mdz2 + mpi2+ \
        2.*np.sqrt(mdz2+vecdz.Mag2())*np.sqrt(mpi2+vecpi.Mag2()) -\
        2.*vecdz.Dot(vecpi)
    return np.sqrt(mdzpi)

#Invariant constants

D0_m = 1864.84 #MeV
triggers = ['D0_Hlt1TrackMVADecision_TOS', 'D0_Hlt1TwoTrackMVADecision_TOS',
            'Hlt2CharmHadDstp2D0Pip_D02KS0KmKp_KS0LLTurboDecision', 'Hlt2CharmHadDstp2D0Pip_D02KS0KmKp_KS0DDTurboDecision',
            'Hlt2CharmHadDstp2D0Pip_D02KS0KmKp_KS0LL_LTUNBTurboDecision', 'Hlt2CharmHadDstp2D0Pip_D02KS0KmKp_KS0DD_LTUNBTurboDecision']

Vars = {'dm':'F',
        'mD0pi':'F',
        'mKSKp':'F',
        'mKSKm':'F',
        'mKK':'F',
        'D0_BPLTIME_norm':'F',
        'Dst_DTF_D0_CTAU_norm':'F',
        'pis_ID':'F',
        'trg':'I',
        'D0_TIP':'F',
        'log_D0_IPCHI2':'F',
        'nPV': 'I'}
# VarsName = Vars.keys() # List of the strings in Vars BUT NOT IN ORDER

t = DefineTree(Vars, 'SkimmedD02KsKK')
iev=0
nevents = 2000
for ev in ch:
    iev+=1
    if iev >nevents: break
    m_pip_pim = GetMass(ev.KS_pip_PX + ev.KS_pim_PX,
                        ev.KS_pip_PY + ev.KS_pim_PY,
                        ev.KS_pip_PZ + ev.KS_pim_PZ,
                        ev.KS_pip_E + ev.KS_pim_E)


    #TRIGGERS

    if (ev.Dst_DTFKS_KS_CTAU > 0 and
        ev.pis_PT > 200 and ev.pis_IPCHI2 < 20 and ev.pis_GHOSTPROB < 0.25 and ev.pis_PIDe < 4 and
        485 < m_pip_pim < 510 and
        -80 < ev.D0_TIP < 80 and #??
        ev.D0_IPCHI2 < 9 and
        1845 < ev.Dst_DTF_D0_M < 1885):

        # #m_D0_pis

        # D0_px = ev.Dst_DTF_D0_KS_PX + ev.Dst_DTF_D0_h1_PX + ev.Dst_DTF_D0_h2_PX
        # D0_py = ev.Dst_DTF_D0_KS_PY + ev.Dst_DTF_D0_h1_PY + ev.Dst_DTF_D0_h2_PY
        # D0_pz = ev.Dst_DTF_D0_KS_PZ + ev.Dst_DTF_D0_h1_PZ + ev.Dst_DTF_D0_h2_PZ
        # D0_E = np.sqrt(D0_m*D0_m + D0_px*D0_px + D0_py*D0_py + D0_pz*D0_pz)
        # m_KS_daup_daum = GetMass(D0_px, D0_py, D0_pz, D0_E)
        #print(m_KS_daup_daum)


        #Calculate the variables that need some work
        m_D0_pis =GetMDzPi(ev, {'Kp':'h2', 'Km': 'h1', 'KS': 'KS'})
        # GetMass(D0_px + ev.Dst_DTF_pis_PX, D0_py + ev.Dst_DTF_pis_PY, D0_pz + ev.Dst_DTF_pis_PZ, D0_E + ev.Dst_DTF_pis_E)

        m_KS_daum = GetMass(ev.Dst_DTF_D0_h1_PX + ev.Dst_DTF_D0_KS_PX,
                            ev.Dst_DTF_D0_h1_PY + ev.Dst_DTF_D0_KS_PY,
                            ev.Dst_DTF_D0_h1_PZ + ev.Dst_DTF_D0_KS_PZ,
                            ev.Dst_DTF_D0_h1_E + ev.Dst_DTF_D0_KS_E)

        m_KS_daup = GetMass(ev.Dst_DTF_D0_h1_PX + ev.Dst_DTF_D0_KS_PX,
                            ev.Dst_DTF_D0_h1_PY + ev.Dst_DTF_D0_KS_PY,
                            ev.Dst_DTF_D0_h1_PZ + ev.Dst_DTF_D0_KS_PZ,
                            ev.Dst_DTF_D0_h1_E + ev.Dst_DTF_D0_KS_E)

        m_daup_daum = GetMass(ev.Dst_DTF_D0_h1_PX + ev.Dst_DTF_D0_h2_PX,
                            ev.Dst_DTF_D0_h1_PY + ev.Dst_DTF_D0_h2_PY,
                            ev.Dst_DTF_D0_h1_PZ + ev.Dst_DTF_D0_h2_PZ,
                            ev.Dst_DTF_D0_h1_E + ev.Dst_DTF_D0_h2_E)
        #Triggers variable
        trg = 0
        for i in range (0, len(triggers)):
            if (ev.__getattr__(triggers[i]) == 1):
                trg += np.power(2, i)

        # t['vars'][VarsName[0]][0] = ev.Dst_DTF_DM
        # t['vars'][VarsName[1]][0] = m_D0_pis
        # t['vars'][VarsName[2]][0] = m_KS_daup
        # t['vars'][VarsName[3]][0] = m_KS_daum
        # t['vars'][VarsName[4]][0] = m_daup_daum
        # t['vars'][VarsName[5]][0] = ev.D0_BPVLTIME / 0.000410
        # t['vars'][VarsName[6]][0] =  # ev.Dst_DTF_D0_CTAU / 0.1229
        # t['vars'][VarsName[7]][0] = ev.pis_ID
        # t['vars'][VarsName[8]][0] = trg
        # t['vars'][VarsName[9]][0] = ev.D0_TIP
        # t['vars'][VarsName[10]][0] = np.log(ev.D0_IPCHI2)

        t['vars']['dm'][0] = ev.Dst_DTF_DM
        # print ('type ev.Dst_DTF_DM = ' , type(ev.Dst_DTF_DM))
        # print ('type m_d0_pis = '  , type(m_D0_pis))
        t['vars']['mD0pi'][0] = m_D0_pis
        t['vars']['mKSKp'][0] = m_KS_daup
        t['vars']['mKSKm'][0] = m_KS_daum
        t['vars']['mKK'][0] = m_daup_daum
        t['vars']['D0_BPLTIME_norm'][0] = ev.D0_BPVLTIME / 0.000410
        t['vars']['Dst_DTF_D0_CTAU_norm'][0] = ev.Dst_DTF_D0_CTAU / 0.1229
        t['vars']['pis_ID'][0] = ev.pis_ID
        t['vars']['trg'][0] = trg
        t['vars']['D0_TIP'][0] = ev.D0_TIP
        t['vars']['log_D0_IPCHI2'][0] = np.log(ev.D0_IPCHI2)
        t['vars']['nPV'][0] = ev.nPV

        t['tree'].Fill()

SaveTree(t['tree'], 'SkimmedD02KsKK_ ' + str(nevents) + '_events.root')

raw_input('Press [Enter] to continue')
