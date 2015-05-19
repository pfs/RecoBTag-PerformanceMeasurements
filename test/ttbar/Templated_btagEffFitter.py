import optparse
import os,sys
import json
import commands
import ROOT

SLICEBINS  = [(20,1000),(20,30),(30,50),(50,80),(80,120),(120,210),(210,320),(320,1000)]
SLICEVAR   = 'jetpt'
SYSTVARS   = ['','jesup','jesdn','jerup','jerdn','puup','pudn']

"""
Project trees from files to build the templates
"""
def prepareTemplates(tagger,taggerDef,var,varRange,inDir,outDir):
    
    print '...starting %s for %s'%(tagger,var)

    #prepare templates
    baseHisto=ROOT.TH1F(var,';Discriminator;Jets',50,varRange[0],varRange[1])
    histos={}
    for flav in ['b','c','other','data']:
        nOPs=len(taggerDef)-2
        for i in xrange(0,nOPs):
            for islice in xrange(0,len(SLICEBINS)):
                for systVar in SYSTVARS:
                    if flav=='data' and var!='' : continue
                    key='%s_pass%d_slice%d%s' % (flav, i, islice,systVar)
                    histos[key]=baseHisto.Clone(key)
                    histos[key].SetDirectory(0)
                    histos[key].Sumw2(0)
    baseHisto.Delete()

    #add files to the corresponding chains
    files = [ f for f in os.listdir(inDir) if '.root' in f ]
    chains={'mc':ROOT.TChain('kin')}
    for f in files: 
        key = 'mc' if 'MC' in f else 'data'
        chains[key].Add(inDir+'/'+f)

    #fill histos
    for key in chains:
        for i in xrange(0,chains[key].GetEntries()):
            chains[key].GetEntry(i)

            for systVar in SYSTVARS:
                
                wgtIdx, systIdx = 0, 0
                if systVar=='jesup' : wgtIdx, systIdx = 1, 1
                if systVar=='jesdn' : wgtIdx, systIdx = 2, 2
                if systVar=='jerup' : wgtIdx, systIdx = 3, 3
                if systVar=='jerdn' : wgtIdx, systIdx = 4, 4
                weight      = chains[key].weight[wgtIdx]

                #no need to proceed if event is not selected
                if weight==0: continue

                #variable to slice on, variable to be fit, and tagger to apply
                sliceVarVal = getattr(chains[key],SLICEVAR)[systIdx] if SLICEVAR=='jetpt' else getattr(chains[key],SLICEVAR) 
                varVal      = getattr(chains[key],var)               if var=='jp'        else getattr(chains[key],var)[systIdx]
                taggerVal   = getattr(chains[key],tagger+'Tagger')

                #determine categories
                passSlice=[]
                for islice in xrange(0,len(SLICEBINS)):
                    if sliceVarVal<=SLICEBINS[islice][0] or sliceVarVal>SLICEBINS[islice][1] : continue
                    passSlice.append(islice)

                #assign flavour
                flav='other'
                if chains[key].flavour==3: flav='b'
                if chains[key].flavour==2: flav='c'
                if key=='data' : flav='data'

                #fill the histos
                for islice in passSlice:
                    hkey='%s_pass0_slice%d%s'%(flav,islice,systVar)
                    histos[hkey].Fill(varVal,weight)
                    for i in xrange(2,len(taggerDef)-1):
                        if taggerVal< taggerDef[i] : continue
                        hkey='%s_pass%d_slice%d%s'%(flav,i-1,islice,systVar)
                        histos[hkey].Fill(varVal,weight)

    #save templates to file
    fOut=ROOT.TFile.Open('%s/%s_templates/%s.root'%(outDir,var,tagger),'RECREATE')
    for key in histos : histos[key].Write()
    fOut.Close()


"""
Wrapper to be used when run in parallel
"""
def runPrepareTemplatesPacked(args):
    tagger, taggerDef, var, varRange, inDir, outDir = args
    try:
        return prepareTemplates(tagger=tagger,
                                taggerDef=taggerDef,
                                var=var,
                                varRange=varRange,
                                inDir=inDir,
                                outDir=outDir)
    except :
        print 50*'<'
        print "  Problem found (%s) baling out of this task" % sys.exc_info()[1]
        print 50*'<'
        return False


"""
steer the script
"""
def main():

    #configuration
    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--taggers',            dest='taggers'  ,          help='json with list of taggers',    default=None,        type='string')
    parser.add_option('-i', '--inDir',              dest='inDir',              help='input directory with files',   default=None,        type='string')
    parser.add_option('-v', '--var',                dest='var',                help='templated variable',           default='kindisc',   type='string')
    parser.add_option(      '--recycleTemplates',   dest='recycleTemplates',   help='recycleTemplates',             default=False,       action='store_true')
    parser.add_option('-n', '--njobs',              dest='njobs',              help='# jobs to run in parallel',    default=0,           type='int')
    parser.add_option('-o', '--outDir',             dest='outDir',             help='output directory',             default='analysis',  type='string')
    (opt, args) = parser.parse_args()
    
    #read list of samples
    taggersFile = open(opt.taggers,'r')
    taggersList=json.load(taggersFile,encoding='utf-8').items()
    taggersFile.close()

    if not opt.recycleTemplates:
        task_list=[]
        for var,varMin,varMax in [('kindisc',-1,1),('jp',0,2),('close_mlj',0,200)]:
            os.system('mkdir -p %s/%s_templates'%(opt.outDir,var))
            for tagger,taggerDef in taggersList:
                if var==tagger : continue
                task_list.append((tagger,taggerDef,var,(varMin,varMax),opt.inDir,opt.outDir))
        #task_list=list(set(task_list))
        print '%s jobs to run in %d parallel threads' % (len(task_list), opt.njobs)
        #run the analysis jobs                                                                                                                                                          
        if opt.njobs == 0:
            for tagger,taggerDef,var,varRange,inDir,outDir in task_list:
                prepareTemplates(tagger=tagger,
                                 taggerDef=taggerDef,
                                 var=var,
                                 varRange=varRange,
                                 inDir=inDir,
                                 outDir=outDir)
        else:
            from multiprocessing import Pool
            pool = Pool(opt.njobs)
            pool.map(runPrepareTemplatesPacked, task_list)

    #all done here
    exit(0)



"""
for execution from another script
"""
if __name__ == "__main__":
    sys.exit(main())
