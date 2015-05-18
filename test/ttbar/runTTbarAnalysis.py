import optparse
import os,sys
import json
import commands
import ROOT
import pickle
import math
from array import array
from storeTools import getEOSlslist

LUMI=100
CHANNELS={-11*11:'ll', -13*13:'ll', -11*13:'emu'}

"""
Perform the analysis on a single file
"""
def runTTbarAnalysis(inFile, outFile, wgt):

    #prepare output and histograms
    outF=ROOT.TFile.Open(outFile,'RECREATE')
    baseHistos={
        'npv'    : ROOT.TH1F('npv',    ';N_{PV,good}-N_{HS};Events',              50, 0, 50),
        'rho'    : ROOT.TH1F('rho',    ';#rho [GeV];Events',                      50, 0, 30),
        'mll'    : ROOT.TH1F('mll',    ';Dilepton invariant mass [GeV];Events',   50, 0, 250),
        'met'    : ROOT.TH1F('met',    ';Missing transverse energy [GeV];Events', 50, 0, 250),
        'jetpt'  : ROOT.TH1F('jetpt',  ';Transverse momentum [GeV]; Events',      50, 0, 200),
        'jeteta' : ROOT.TH1F('jeteta', ';Pseudo-rapidity; Events',                25, 0, 2.5),
        'njets'  : ROOT.TH1F('njets',  ';Jet multiplicity;Events;',               6,  2, 8),
        'jp'     : ROOT.TH1F('jp',     ';Jet probability; Jets;',                 50, 0, 2.5)
        }

    kinTree=ROOT.TTree('kin','kin training')
    kinTree.SetDirectory(outF)
    kinVars={'signal'     : (array( 'i', [ 0 ]), ';b-jet;Jets',                         2,  0, 2 ),
             'close_mlj'  : (array( 'f', [ 0.]), ';M(lepton,jet) [GeV]; Jets',          50, 0, 250),
             'close_deta' : (array( 'f', [ 0.]), ';#Delta#eta(lepton,jet); Jets',       50, 0, 4),
             'close_dphi' : (array( 'f', [ 0.]), ';#Delta#phi(lepton,jet) [rad]; Jets', 50, 0, 3.15),
             'close_ptrel': (array( 'f', [ 0.]), ';p_{T}^{rel}(lepton,jet) [GeV];Jets', 50, 0, 1),
             'far_mlj'    : (array( 'f', [ 0.]), ';M(lepton,jet) [GeV]; Jets',          50, 0, 250),
             'far_deta'   : (array( 'f', [ 0.]), ';#Delta#eta(lepton,jet); Jets',       50, 0, 4),
             'far_dphi'   : (array( 'f', [ 0.]), ';#Delta#phi(lepton,jet) [rad]; Jets', 50, 0, 3.15),
             'far_ptrel'  : (array( 'f', [ 0.]), ';p_{T}^{rel}(lepton,jet) [GeV];Jets', 50, 0, 1) }
    for v in kinVars:
        vtype = 'I' if v=='signal' else 'F'
        kinTree.Branch( v, kinVars[v][0], '%s/%s' % (v,vtype) )
        baseHistos[v] = ROOT.TH1F(v,kinVars[v][1],kinVars[v][2],kinVars[v][3],kinVars[v][4])

    #replicate histos per channel
    histos={}
    for ch in CHANNELS.itervalues():
        for key in baseHistos:
            tag='%s_%s' % (ch,key)
            histos[tag]=baseHistos[key].Clone(tag)
            histos[tag].Sumw2()
            histos[tag].SetDirectory(outF)

    #loop over events
    inF=ROOT.TFile.Open(inFile)
    tree=inF.Get("btagana/ttree")
    nentries=tree.GetEntriesFast()
    print '....opened %s -> analysing %d events -> %s' % (inFile,nentries,outFile)
    for i in xrange(0,nentries): 
        tree.GetEntry(i)

        #channel name
        ch=''
        try:
            ch=CHANNELS[tree.ttbar_chan]
        except:
            continue

        if tree.ttbar_nl<2 : continue
        lp4=[]
        for il in xrange(0,tree.ttbar_nl):
            lp4.append( ROOT.TLorentzVector() )
            lp4[-1].SetPtEtaPhiM(tree.ttbar_lpt[il],tree.ttbar_leta[il],tree.ttbar_lphi[il],0.)
        mll=(lp4[0]+lp4[1]).M()


        selJets=[]
        ljkinVars=[]
        for ij in xrange(0,tree.nJet):
            if tree.Jet_pt[ij]<30 or math.fabs(tree.Jet_eta[ij])>2.5 : continue
            jp4=ROOT.TLorentzVector()
            jp4.SetPtEtaPhiM(tree.Jet_pt[ij],tree.Jet_eta[ij],tree.Jet_phi[ij],0.)

            minDRlj=9999.
            ljkin=[]
            for il in xrange(0,tree.ttbar_nl):
                dphi=ROOT.TMath.Abs(lp4[il].DeltaPhi(jp4))
                deta=ROOT.TMath.Abs(lp4[il].Eta()-jp4.Eta())
                dr=lp4[il].DeltaR(jp4)
                minDRlj=ROOT.TMath.Min( minDRlj, dr )
                ptrel=ROOT.Math.VectorUtil.Perp(lp4[il].Vect(),jp4.Vect().Unit())/lp4[il].P();
                ljkin.append( ( (lp4[il]+jp4).M(),dr,deta,dphi,ptrel ) )                
            if minDRlj<0.4 : continue

            selJets.append(ij)
            ljkin=sorted(ljkin, key=lambda pair: pair[1])
            ljkinVars.append( ljkin )

        #require two jets
        njets=len(selJets)
        if njets<2 : continue
        histos[ch+'_njets'].Fill(njets,wgt)

        #n-1 plots
        zCand   = True if 'll' in ch and ROOT.TMath.Abs(mll-91)<15 else False
        passMet = True if 'emu' in ch or tree.ttbar_metpt>40 else False
        if passMet   : histos[ch+'_mll'].Fill(mll,wgt)
        if not zCand : histos[ch+'_met'].Fill(tree.ttbar_metpt,wgt)
        if zCand : continue
        if not passMet : continue
        histos[ch+'_npv'].Fill(tree.nPV-1,wgt)
        histos[ch+'_rho'].Fill(tree.ttbar_rho,wgt)        

        #jet control plots
        for ij in selJets:
            histos[ch+'_jetpt'].Fill(tree.Jet_pt[ij],wgt)
            histos[ch+'_jeteta'].Fill(ROOT.TMath.Abs(tree.Jet_eta[ij]),wgt)
            histos[ch+'_jp'].Fill(ROOT.TMath.Abs(tree.Jet_Proba[ij]),wgt)
            kinVars['signal'][0][0]     = 1 if ROOT.TMath.Abs(tree.Jet_flavour[ij])==5 else 0
            kinVars['close_mlj'][0][0]  = ljkinVars[ij][0][0]
            kinVars['close_deta'][0][0] = ljkinVars[ij][0][2]
            kinVars['close_dphi'][0][0] = ljkinVars[ij][0][3]
            kinVars['close_ptrel'][0][0]= ljkinVars[ij][0][4]
            kinVars['far_mlj'][0][0]    = ljkinVars[ij][1][0]
            kinVars['far_deta'][0][0]   = ljkinVars[ij][1][2]
            kinVars['far_dphi'][0][0]   = ljkinVars[ij][1][3]
            kinVars['far_ptrel'][0][0]  = ljkinVars[ij][1][4]
            for v in kinVars : histos[ch+'_'+v].Fill(kinVars[v][0][0],wgt)
            kinTree.Fill()

    inF.Close()
    
    #dump results to file
    outF.cd()
    for h in histos: histos[h].Write()
    kinTree.Write()                            
    outF.Close()


"""
Wrapper to be used when run in parallel
"""
def runTTbarAnalysisPacked(args):
    inFile, outFile, wgt = args
    try:
        return runTTbarAnalysis(inFile=inFile, outFile=outFile, wgt=wgt)
    except :
        print 50*'<'
        print "  Problem  (%s) with %s continuing without"%(sys.exc_info()[1],inFile)
        print 50*'<'
        return False


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',      default=None,        type='string')
    parser.add_option('-i', '--inDir',       dest='inDir',       help='input directory with files',   default=None,        type='string')
    parser.add_option('-l', '--lumi',        dest='lumi',        help='integrated luminosity to use', default=100.,        type='float')
    parser.add_option('-o', '--outDir',      dest='outDir',      help='output directory',             default='analysis',  type='string')
    parser.add_option('-n', '--njobs',       dest='njobs',       help='# jobs to run in parallel',    default=0,           type='int')
    (opt, args) = parser.parse_args()

    #update luminosity
    global LUMI
    LUMI=opt.lumi
    
    #read list of samples
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile,encoding='utf-8').items()
    jsonFile.close()

    #prepare output
    if len(opt.outDir) : os.system('mkdir -p %s' % opt.outDir)

    #read normalization
    xsecWgts, integLumi = {}, {}
    cache='%s/.xsecweights.pck'%opt.outDir
    try:
        cachefile = open(cache, 'r')
        xsecWgts  = pickle.load(cachefile)
        integLumi = pickle.load(cachefile)
        cachefile.close()        
        print 'Normalization read from cache (%s)' % cache
    except:
        print 'Computing original number of events and storing in cache, this may take a while if it\'s the first time'
        for tag,sample in samplesList: 

            if sample[1]==1 : 
                xsecWgts[tag]=1.0
                continue

            input_list=getEOSlslist(directory=opt.inDir+'/'+tag)            
            xsec=sample[0]            
            norigEvents=0
            for f in input_list:
                fIn=ROOT.TFile.Open(f)
                norigEvents+=fIn.Get('allEvents/hEventCount').GetBinContent(1)
                fIn.Close()
            xsecWgts[tag]  = LUMI*xsec/norigEvents if norigEvents>0 else 0
            integLumi[tag] = norigEvents/xsec      if norigEvents>0 else 0
            print '... %s cross section=%f pb #orig events=%d lumi=%3.2f/fb' % (tag,xsec,norigEvents,integLumi[tag]/1000.)

        #dump to file
        cachefile=open(cache,'w')
        pickle.dump(xsecWgts, cachefile, pickle.HIGHEST_PROTOCOL)
        pickle.dump(integLumi, cachefile, pickle.HIGHEST_PROTOCOL)
        cachefile.close()
        print 'Produced normalization cache (%s)'%cache

    #create the analysis jobs
    task_list = []
    for tag,_ in samplesList:
        input_list=getEOSlslist(directory=opt.inDir+'/'+tag)
        wgt = xsecWgts[tag]
        for nf in xrange(0,len(input_list)) : 
            outF='%s/%s_%d.root'%(opt.outDir,tag,nf)
            task_list.append( (input_list[nf],outF,wgt) )
    task_list=list(set(task_list))
    print '%s jobs to run in %d parallel threads' % (len(task_list), opt.njobs)

    #run the analysis jobs
    if opt.njobs == 0:
        for inFile, outFile,wgt in task_list: 
            runTTbarAnalysis(inFile=inFile, outFile=outFile, wgt=wgt)
    else:
        from multiprocessing import Pool
        pool = Pool(opt.njobs)
        pool.map(runTTbarAnalysisPacked, task_list)

    #merge the outputs
    for tag,_ in samplesList:
        os.system('hadd -f %s/%s.root %s/%s_*.root' % (opt.outDir,tag,opt.outDir,tag) )
        os.system('rm %s/%s_*.root' % (opt.outDir,tag) )
    print 'Analysis results are available in %s' % opt.outDir

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
