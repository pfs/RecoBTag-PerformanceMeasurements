import optparse
import os,sys
import json
import commands
import ROOT
import pickle
import math

LUMI=100

from storeTools import getEOSlslist

"""
Perform the analysis on a single file
"""
def runTTbarAnalysis(inFile, outFile, wgt):

    #prepare output
    outF=ROOT.TFile.Open(outFile,'RECREATE')
    histos={
        'njets' : ROOT.TH1F('njets', ';Jet multiplicity;Events;', 6,  0, 6),
        'jp'    : ROOT.TH1F('jp',    ';Jet probability; Jets;',   50, 0, 2.5)
        }
    for h in histos:
        histos[h].Sumw2()
        histos[h].SetDirectory(outF)
    
    #loop over events
    inF=ROOT.TFile.Open(inFile)
    tree=inF.Get("btagana/ttree")
    nentries=tree.GetEntriesFast()
    print '....opened %s -> analysing %d events -> %s' % (inFile,nentries,outFile)
    for i in xrange(0,nentries): 
        tree.GetEntry(i)
        nselJets=0
        for ij in xrange(0,tree.nJet):
            if tree.Jet_pt[ij]<30 or math.fabs(tree.Jet_eta[ij])>2.4 : continue
            jp4=ROOT.TLorentzVector()
            jp4.SetPtEtaPhiM(tree.Jet_pt[ij],tree.Jet_eta[ij],tree.Jet_phi[ij],0.)

            minDRlj=9999.
            for il in xrange(0,tree.ttbar_nl):
                lp4=ROOT.TLorentzVector()
                lp4.SetPtEtaPhiM(tree.ttbar_lpt[il],tree.ttbar_leta[il],tree.ttbar_lphi[il],0.)
                minDRlj=ROOT.TMath.Min(minDRlj,lp4.DeltaR(jp4))
            if minDRlj<0.4 : continue

            nselJets+=1
            histos['jp'].Fill(tree.Jet_Proba[ij],wgt)

        histos['njets'].Fill(nselJets,wgt)

    inF.Close()
    
    #dump results to file
    outF.cd()
    for h in histos: histos[h].Write()
    outF.Close()


"""
Wrapper to be used when run in parallel
"""
def runTTbarAnalysisPacked(args):
    inFile, outFile, wgt = args
    try:
        return runTTbarAnalysis(inFile=inFile, outFile=outFile, wgt=wgt)
    except:
        print 50*'<'
        print "  Problem with %s continuing without"%inFile
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
    xsecWgts={}
    cache='%s/.xsecweights.pck'%opt.outDir
    try:
        cachefile = open(cache, 'r')
        xsecWgts = pickle.load(cachefile)
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
            print '... %s cross section=%f pb #orig events=%d ' % (tag,xsec,norigEvents)
            xsecWgts[tag]=xsec/norigEvents if norigEvents>0 else 0
        
        #dump to file
        cachefile=open(cache,'w')
        pickle.dump(xsecWgts, cachefile, pickle.HIGHEST_PROTOCOL)
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
