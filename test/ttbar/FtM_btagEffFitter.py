import optparse
import os,sys
import json
import commands
import ROOT
import pickle
from plotter import Plot

RELTAGEFFVAR  = 0.5
JETMULTCATEGS = [2,3,4]
SLICEBINS     = [(30,60),(60,120),(120,320)]
SLICEVAR      = 'jetpt'
SYSTVARS      = ['nom',
                 'jesup','jesdn',
                 'jerup','jerdn',
                 'puup','pudn',
                 'trigup','trigdn',
                 'selup','seldn',
                 'qcdscaledn','qcdscaleup',
                 'hdampdn','hdampup'                
                 ]

"""
Project trees from files to build the templates
"""
def prepareTemplates(sampleTag,tagger,taggerDef,sliceCat,channelList,btagSFTool,expEffUrl,inDir,outDir):

    baseName='%s_count_%s_%d%d'%(tagger,SLICEVAR,sliceCat[0]++1,sliceCat[1]+1)

    #prepare task systematics and tag counting
    nOPs=len(taggerDef)-2
    tagCounting={'nom':[0]*(nOPs-1)}
    TASKSYSTVARS=SYSTVARS[:]
    effCategs=[0]
    effCategs.append(sliceCat[0]+1)
    if sliceCat[0]!=sliceCat[1]:
        effCategs.append(sliceCat[1]+1)
    for idx in effCategs:
        for eff in ['beff','ceff','leff']:
            for sign in ['Up','Down']:
                TASKSYSTVARS.append( '%s%d%s'%(eff,idx,sign) )
                tagCounting['%s%d%s'%(eff,idx,sign)]=[0]*(nOPs-1)
        
    #prepare histograms
    histos={}
    nJetMultCategs=len(JETMULTCATEGS)
    histoVars=TASKSYSTVARS
    for iop in xrange(1,nOPs):
        for var in histoVars:
            name='%s_pass%d_%s'%(baseName,iop,var)
            histos[name]=ROOT.TH1F(name,';Category;Events',3*nJetMultCategs,0,3*nJetMultCategs)
            for ij in xrange(0,nJetMultCategs):
                histos[name].GetXaxis().SetBinLabel(ij*3+1,'%dj0t'%JETMULTCATEGS[ij])
                histos[name].GetXaxis().SetBinLabel(ij*3+2,'%dj1t'%JETMULTCATEGS[ij])
                histos[name].GetXaxis().SetBinLabel(ij*3+3,'%dj2t'%JETMULTCATEGS[ij])
            histos[name].Sumw2()
            histos[name].SetDirectory(0)

    #report
    print '...starting %s for category %3.0f to %3.0f / %3.0f to %3.0f'%(tagger,
                                                                         SLICEBINS[sliceCat[0]][0],SLICEBINS[sliceCat[0]][1],
                                                                         SLICEBINS[sliceCat[1]][0],SLICEBINS[sliceCat[1]][1])
    print '   base name is %s'%baseName
    print '   task systematics are',TASKSYSTVARS
    print '   will fill %d histograms'%len(histos)
    print '   for events in',channelList

    #read expected efficiencies from cache
    cachefile=open(expEffUrl,'r')
    expEffGr=pickle.load(cachefile)
    cachefile.close()

    #fill histos
    chain=ROOT.TChain('ftm')
    chain.Add('%s/%s.root'%(inDir,sampleTag))
    totalEntries=chain.GetEntries()
    for i in xrange(0,totalEntries):

        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )
            
        chain.GetEntry(i)

        #require matching channel
        if not chain.ttbar_chan in channelList : continue

        #require jet multiplicity within the ones selected
        jetMultCateg=-1
        for ijm in xrange(0,len(JETMULTCATEGS)):
            if  chain.jetmult != JETMULTCATEGS[ijm] : continue
            jetMultCateg=ijm
            break
        if jetMultCateg<0: continue
            
        #check if jets are in required category
        sliceVarVal1=getattr(chain,SLICEVAR)[0]
        sliceVarVal2=getattr(chain,SLICEVAR)[1]
        jetCats=[-1,-1]
        if sliceVarVal1>SLICEBINS[sliceCat[0]][0] and sliceVarVal1<SLICEBINS[sliceCat[0]][1]:
            jetCats[0]=sliceCat[0]+1
            if sliceVarVal2>SLICEBINS[sliceCat[1]][0] and sliceVarVal2<SLICEBINS[sliceCat[1]][1]:
                jetCats[1]=sliceCat[1]+1
        if sliceVarVal1>SLICEBINS[sliceCat[1]][0] and sliceVarVal1<SLICEBINS[sliceCat[1]][1]:
            jetCats[0]=sliceCat[1]+1
            if sliceVarVal2>SLICEBINS[sliceCat[0]][0] and sliceVarVal2<SLICEBINS[sliceCat[0]][1]:
                jetCats[1]=sliceCat[0]+1
        if jetCats[0]<0 or jetCats[1]<0: continue
                
        #reset tag counters
        for tagVar in tagCounting: tagCounting[tagVar]=[0]*(nOPs-1)
        for ij in xrange(0,2):
                    
            #efficiency variation categories
            effCategs=[0]
            effCategs.append( jetCats[ij])

            #jet pT
            jpt=chain.jetpt[ij]

            #MC truth
            flav,flavName=abs(chain.flavour[ij]),'udsg'
            if flav==5: flavName='b'
            if flav==4: flavName='c'

            #tagger value
            taggerVal = getattr(chain,tagger)[ij]

            #count tags
            for iop in xrange(1,nOPs):

                expEff = 1.0
                try:
                    expEff=expEffGr[tagger][flavName][iop].Eval(jpt) 
                except:
                    pass
                
                passTag = False if taggerVal<taggerDef[iop+1] else True
                tagCounting['nom' ][iop-1] += int(passTag)
                    
                passTagBeffUp,passTagBeffDown = passTag,passTag
                if flav==5:
                    passTagBeffUp=btagSFTool.modifyBTagsWithSF(passTagBeffUp,1.0*(1+RELTAGEFFVAR),expEff)
                    passTagBeffDown=btagSFTool.modifyBTagsWithSF(passTagBeffDown,1.0*(1-RELTAGEFFVAR),expEff)
                for idx in effCategs:
                    tagCounting['beff%dUp'%idx][iop-1] += int(passTagBeffUp)
                    tagCounting['beff%dDown'%idx][iop-1] += int(passTagBeffDown)


                passTagCeffUp,passTagCeffDown = passTag,passTag
                if flav==4:
                    passTagCeffUp=btagSFTool.modifyBTagsWithSF(passTagCeffUp,1.0*(1+RELTAGEFFVAR),expEff)
                    passTagCeffDown=btagSFTool.modifyBTagsWithSF(passTagCeffDown,1.0*(1-RELTAGEFFVAR),expEff)
                for idx in effCategs:
                    tagCounting['ceff%dUp'%idx][iop-1] += int(passTagCeffUp)
                    tagCounting['ceff%dDown'%idx][iop-1] += int(passTagCeffDown)

                passTagLeffUp,passTagLeffDown = passTag,passTag
                if flav!=4 and flav!=5:
                    passTagLeffUp=btagSFTool.modifyBTagsWithSF(passTagLeffUp,1.0*(1+RELTAGEFFVAR),expEff)
                    passTagLeffDown=btagSFTool.modifyBTagsWithSF(passTagLeffDown,1.0*(1-RELTAGEFFVAR),expEff)
                for idx in effCategs:
                    tagCounting['leff%dUp'%idx][iop-1] += int(passTagLeffUp)
                    tagCounting['leff%dDown'%idx][iop-1] += int(passTagLeffDown)

        #print sliceVarVal1,sliceVarVal2,'->',effCategs
        #print 'Flavours:',chain.flavour[0],chain.flavour[1]
        #print tagCounting

        #fill the histograms        
        for iop in xrange(1,nOPs):
            for var in histoVars:
                name='%s_pass%d_%s'%(baseName,iop,var)
                
                weight=chain.weight[0]
                if var=='jesup':      weight=chain.weight[1]
                if var=='jesdn':      weight=chain.weight[2]
                if var=='jerup':      weight=chain.weight[3]
                if var=='jerdn':      weight=chain.weight[4]
                if var=='puup':       weight=chain.weight[5]
                if var=='pudn':       weight=chain.weight[6]
                if var=='trigup':     weight=chain.weight[7]
                if var=='trigdn':     weight=chain.weight[8]
                if var=='selup':      weight=chain.weight[9]
                if var=='seldn':      weight=chain.weight[10]
                if var=='qcdscaleup': weight=chain.weight[11]
                if var=='qcdscaledn': weight=chain.weight[12]
                if var=='hdampup':    weight=chain.weight[13]
                if var=='hdampdn':    weight=chain.weight[14]
                nbtags=tagCounting['nom'][iop-1]
                if 'beff' in var or 'ceff' in var or 'leff' in var:
                    nbtags=tagCounting[var][iop-1]
                histos[name].Fill(nbtags+jetMultCateg*3,weight)
                #print name,nbtags,weight
        #print tagCounting
        #raw_input()
                     
    #save templates to file
    fOut=ROOT.TFile.Open('%s/FtM/%s.root'%(outDir,sampleTag),'UPDATE')
    print '    output stored in %s'%fOut.GetName()
    for key in histos : histos[key].Write()
    fOut.Close()



"""
Wrapper to be used when run in parallel
"""
def runPrepareTemplatesPacked(args):
    sampleTag,tagger,taggerDef,sc,channelList,btagSFTool,inDir,outDir,expEffUrl in args
    try:
        return prepareTemplates(sampleTag=sampleTag,
                                tagger=tagger,
                                taggerDef=taggerDef,
                                sliceCat=sc,
                                channelList=channelList,
                                btagSFTool=btagSFTool,
                                expEffUrl=expEffUrl,
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
    parser.add_option('-j', '--json',        dest='json'  ,      help='json with list of files',      default=None,        type='string')
    parser.add_option('-i', '--inDir',              dest='inDir',              help='input directory with files',   default=None,        type='string')
    parser.add_option('-e',  '--exp',               dest='expEff',             help='pickle file with the expected efficiencies',   default='data/.expEff.pck',      type='string')
    parser.add_option('-n', '--njobs',              dest='njobs',              help='# jobs to run in parallel',    default=0,           type='int')
    parser.add_option('-o', '--outDir',             dest='outDir',             help='output directory',             default='analysis',  type='string')
    parser.add_option(      '--channels',           dest='channels',           help='channels to use',              default='-121,-143,-169',  type='string')
    parser.add_option(      '--recycleTemplates',   dest='recycleTemplates',   help='do not regenerate templates',  default=False,       action='store_true')
    (opt, args) = parser.parse_args()

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gSystem.CompileMacro("BtagUncertaintyComputer.cc","fkgd","libBtagUncertaintyComputer");
    ROOT.gSystem.Load("libBtagUncertaintyComputer.so")
    btagSFTool=ROOT.BTagSFUtil()

    #channels to filter
    channelList=[ int(k) for k in opt.channels.split(',') ]

    #read list of samples
    jsonFile = open(opt.json,'r')
    samplesList=json.load(jsonFile,encoding='utf-8').items()
    jsonFile.close()
    
    #read list of taggers
    taggersFile = open(opt.taggers,'r')
    taggersList=json.load(taggersFile,encoding='utf-8').items()
    taggersFile.close()
   
    #prepare output
    os.system('mkdir -p %s/FtM'%(opt.outDir))

    #all possible slice categs
    sliceCategs=[]
    for i in xrange(0,len(SLICEBINS)):
        for j in xrange (i,len(SLICEBINS)):
            sliceCategs.append((i,j))

    #re-create templates 
    if not opt.recycleTemplates:
   
        task_list=[]
        for sampleTag,_ in samplesList:

            #preapare output file
            fOut=ROOT.TFile.Open('%s/FtM/%s.root'%(opt.outDir,sampleTag),'RECREATE')
            fOut.Close()

            for tagger,taggerDef in taggersList:
                for sc in sliceCategs:
                    task_list.append((sampleTag,tagger,taggerDef,sc,channelList,btagSFTool,opt.inDir,opt.outDir,opt.expEff))
    
        print '%s jobs to run in %d parallel threads' % (len(task_list), opt.njobs)
        if opt.njobs == 0:
            for sampleTag,tagger,taggerDef,sc,channelList,btagSFTool,inDir,outDir,expEffUrl in task_list:
                prepareTemplates(sampleTag=sampleTag,
                                 tagger=tagger,
                                 taggerDef=taggerDef,
                                 sliceCat=sc,
                                 channelList=channelList,
                                 btagSFTool=btagSFTool,
                                 expEffUrl=expEffUrl,
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
