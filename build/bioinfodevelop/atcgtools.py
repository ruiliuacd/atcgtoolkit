#!/usr/bin/env python3
'''
Created on 2024年8月24日

@author: pc
'''
import concurrent.futures

import  os,re,argparse,pickle#,time
import config
from bioinfodevelop.analysisUtils.genomeUtils import store_sequence,count_gap
# from bioinfodevelop.analysisUtils.variantUtils import VCF_Data
from bioinfodevelop.analysisUtils import variantUtils
from config import random_str

print("""
*******************************************************************
* atcgtools
* version v0.1 Linux
* Built at Nov 15 2024 21:14:25, by pysam0.19, py3.6
* (C) 2024-present, Liu Lab, x University
* Please report bugs to Rui Liu <ruiliugenetic@nwu.edu.cn>
*******************************************************************

""")

if __name__ == '__main__':
    # path=sys.argv[1] if len(sys.argv)>=2 else os.getcwd()
    parser = argparse.ArgumentParser(description="A multi-command tool")
    subparsers = parser.add_subparsers(dest="analysistype")
    
    countgaps_parser = subparsers.add_parser("countgaps", help="countgaps")
    
    inittopleveltable_parser=subparsers.add_parser('inittopleveltable',help='step 1')
    
    fillcontextOutgroup_parser=subparsers.add_parser('fillcontextNoutgroup',help=""" step 2
        -f exist -> fill context by extract flank seqs from ref.fa
        -a exist -> add column in toplevel and fill outgroup's base to each rec of toplevel table; not increase recs (i.e. SNPs)
        both exist, do the both
    """)
    fillcontextOutgroup_parser.add_argument("-C", "--chromlistfilename", dest="chromlistfilename", help="bed like file , or length file e.g.line: chr1 3424123")
    fillcontextOutgroup_parser.add_argument("-r", "--rflankcontext",nargs=2, help="'Reference genome file' 'flanklen'; recomend 'flanklent'=3, if 'flanklen'>7, file context of db.topvtable and produce fa file store flank seqs ")
    fillcontextOutgroup_parser.add_argument("-a", "--AncensAllelfromVCFBAM",action="append", help="add ancestral allele according outgroup vcf, serve as archic pop. using .VCFBAMconfig file")
    
    vjointT_parser = subparsers.add_parser("VJtransf", help="joint multiple(or just one) vcf files and transform into other variant file type")
    vjointT_parser.add_argument("-v","--vcffile",action="append",help="GATK or samtools/bcftools produced")
    vjointT_parser.add_argument("-o","--outputpre",help="")
    vjointT_parser.add_argument("-c","--chrommap",help="change chr name correspondingly")
    vjointT_parser.add_argument("-d","--dilutetodensity",default=None,help="")
    vjointT_parser.add_argument("-t","--nt",default=1,help="number of chromosomes processed for one time, simultaneously")
    vjointT_parser.add_argument("-r","--r2",default=None,help="prune SNP according r2 invoke plink")
    vjointT_parser.add_argument("-f","--outfmt",help="pedmap/genosnpind/[vcf] indicating that output with .map .ped/.geno .snp .ind/(.vcf unavailable so for) format")
    
    vjointGenedadi_parser = subparsers.add_parser("VJgendadi", help="Generate dadi ds infile")
    vjointGenedadi_parser.add_argument("-q","--quantizingvcf",dest="quantizingvcf",action="append",nargs=2,help="select from some first param of -v ")
    vjointGenedadi_parser.add_argument("-n","--noofindvds2quantizing",dest="noofindvds2quantizing",nargs=2,help="together with -q; first value: noofindvds2quantizing,second value: indicate whether threshold for each pool vcf is important [critical/optional] ")
    vjointGenedadi_parser.add_argument("-R","--reffa",dest="reffa")
    vjointGenedadi_parser.add_argument("-c","--chromlist",dest="chromlist")
    vjointGenedadi_parser.add_argument("-d","--snpperkb",dest="snpperkb")
    vjointGenedadi_parser.add_argument("-i","--includeSNPID",dest="includeSNPID",help="conflict with -d, extract those snpID from -c regions")
    vjointGenedadi_parser.add_argument("-v", "--vcffile", dest="vcffile",action="append",default=[],nargs=2,help="vcffile minAN for indvd; vcffile treatAN for pool. minAN used for filter AN, treatAN used as AN and filter (*1.5<)DP")#indvdNo_Of_POOL=10
    vjointGenedadi_parser.add_argument("-o","--outputfilename",dest="outputfilename")
    vjointGenedadi_parser.add_argument("-a","--alterpresenceForspecifichr",dest="alterpresence",default=[],action="append",help="default not update flags at all(fastest). single 'all':update flags for all exist SNPs in DB; or single 'updatedNnewSNPonly' only update flags for those records execute in insertorUpdatetopleveltable(update or insert SNPs),or a list of -a indicate updating which chr or any str not the chr name indicating don't updat flag")
    vjointGenedadi_parser.add_argument("-t","--nt",default=1,help="number of chromosomes processed for one time, simultaneously")
    
    mutect2_parser = subparsers.add_parser("test", help="Run tutect2")
    mutect2_parser.add_argument("-R", "--reference", required=True, help="Reference genome file")
    
    # args = parser.parse_args()
    (options, args) = parser.parse_known_args()

    print(options, args)
    if options.analysistype.lower()=="countgaps":
        fafname=args[0] if len(args)==1 else 'ragtag.scaffold.fasta'
        if not os.path.exists(fafname):
            fafname=os.path.join(os.getcwd(),fafname)
        print(fafname)
        sequence=store_sequence(fafname)
        # count_gap.py 作用是统计这一步产生的各scaffold上的gap数量
        
        gap = []    # 用于存储序列gap信息的字典    
        with concurrent.futures.ProcessPoolExecutor(max_workers=36) as executor:
            gap = list(executor.map(count_gap,  [sequence[header] for header in sorted(sequence.keys())]) )

        idx=0
        for header in sorted(sequence.keys()):
            counts = gap[idx]
            print(f'scaffold：{header}\tgap数量:{counts}')
            idx+=1
    
    elif options.analysistype=="inittopleveltable":
        rightmostpoponly=True;vcffile=args[0]#sys.argv[2]
        vattools=variantUtils.AncestralAlleletabletools()
        vattools.createtopleveltable(config.topleveltableofvdb,vcffile,rightmostpoponly)
    
    elif options.analysistype=="fillcontextNoutgroup":
        import pysam
        ancestralalleletabletools=variantUtils.AncestralAlleletabletools(database=config.variantsdbname, ip=config.ip)
        outfile=None;flanklen=0;duckrefhandle=None;duckrefindex=None
        if options.rflankcontext:
            flanklen=int(options.rflankcontext[1])
            print(options.rflankcontext,"\noptions:",options)
            if flanklen>7:
                outfile=open(options.rflankcontext[0]+"_snpflankseq.fa",'w')
            duckrefhandle=open(options.rflankcontext[0],'r')
            try:
                duckrefindex = pickle.load(open(options.rflankcontext[0] + ".myfasteridx", 'rb'))
        #             originalspeciesindex = pickle.load(open(originalspeciesref + ".myindex", 'rb'))
            except IOError:
                from bioinfodevelop.analysisUtils import Utils 
                Utils.generateFasterRefIndex(options.rflankcontext[0], options.rflankcontext[0] + ".myfasteridx")
                duckrefindex = pickle.load(open(options.rflankcontext[0] + ".myfasteridx", 'rb'))   
        if options.AncensAllelfromVCFBAM:
            vcfnameKEY_vcfobj_pyBAMfilesVALUE={}
            print('read and rewrited (append) ancestralVcfConfig file. add col when it dose exist, fill alleles from each of',options.AncensAllelfromVCFBAM)
            
            for ancVc in options.AncensAllelfromVCFBAM:
                vcfBconfig=open(ancVc,"r")
                for line in vcfBconfig:
                    vcffilename_obj=re.search(r"vcffilename=(.*)",line.strip())
                    if vcffilename_obj!=None:
                        vcfname=vcffilename_obj.group(1).strip()
                        vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname]=[]
                        vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(variantUtils.VCF_Data(vcfname))
                    elif line.split():
                        vcfnameKEY_vcfobj_pyBAMfilesVALUE[vcfname].append(pysam.AlignmentFile(line.strip(),'rb'))
                vcfBconfig.close()
                # Add a new section or option if it doesn't exist
                # if not config.cfparser.has_section("newsection"):
                #     config.cfparser.add_section("newsection")
                # cfparser.set("newsection", "newoption", "newvalue")
                new_value = config.outgroupVCFBAMconfig_ZJU1ref + ";"+ancVc
                if ancVc not in config.outgroupVCFBAMconfig_ZJU1ref:config.cfparser.set("variantdata", "outgroupVCFBAMconfig_ZJU1ref", new_value.strip(";"))
                # Write changes back to the config file
            with open(os.path.join(config.prjpath,"BDRobot/com/config.properties"), "w") as configfile:
                config.cfparser.write(configfile)
                
        print("fill ancestrallele base type or flank context or both",options.chromlistfilename)
        chromlist=[];chrom_lenlist=[]
        chromlistfile=open(options.chromlistfilename,"r")
        for chrrow in chromlistfile:
            chrrowlist=re.split(r'\s+',chrrow.strip())
            chromlist.append(chrrowlist[0].strip())
            if len(chrrowlist)>=3 and options.chromlistfilename.endswith(".bed"):
                chrom_lenlist.append((chrrowlist[0].strip(),int(chrrowlist[2].strip())))
            else:
                chrom_lenlist.append((chrrowlist[0].strip(),int(chrrowlist[1].strip())))
        ancestralalleletabletools.fillAncestral_context(vcfnameKEY_vcfobj_pyBAMfilesVALUE, chrom_lenlist, duckrefhandle, duckrefindex, outfile, flanklen, config.topleveltableofvdb)
        if options.rflankcontext:
            duckrefhandle.close()
            if flanklen>7:
                outfile.close()
        chromlistfile.close();print("done")
    elif options.analysistype=="VJtransf":#directly copy from basicutil/src/Vcf2Ped.py
  build   from functools import reduce
        from multiprocessing import Process, Manager,Queue
        
        chromchangemap={}
        if options.chrommap:
            f=open(options.chrommap,'r')
            for line in f:
                linelist=re.split(r'\s+',line.strip())
                chromchangemap[linelist[0].strip()]=linelist[1].strip()
            f.close()
        winSNPs=100;stepSNPs=20
        def Vcf2Pedrandomdilut_Wapper(shared_dict,k,rest_queue,rest_listError):
            # import struct
            
            vobj_l=shared_dict[k]['vcfcontentgen']
            tlst=shared_dict[k]['title'][:9]
            if len(vobj_l)>1:
                vcf_aligned="aligned";listOfpopvcfRecsMapByChr=[]
                for vobj in vobj_l[:]:
                    listOfpopvcfRecsMapByChr.append({k:vobj.getVcfListByChrom(k,MQfilter=None)})
                    tlst+=vobj.VcfIndexMap['title'][9:]
                fulloutjoinSNPs=variantUtils.alignmultPopSnpPos(listOfpopvcfRecsMapByChr, "o",None,False)
                shared_dict[k]['vcfcontentgen']=fulloutjoinSNPs[k]
                shared_dict[k]['title']=tlst
            else:
                vcf_aligned="vcf"
                shared_dict[k]['vcfcontentgen']=vobj_l[0].getVcfListByChrom(k,MQfilter=None)
                shared_dict[k]['title']=vobj_l[0].VcfIndexMap['title']
            allchromOrder=shared_dict[k]['chrom']
            print("included all samples:",shared_dict[k]['title'])# print(len(shared_dict[k]['vcfcontent']),type(shared_dict[k]['vcfcontent']))
            shared_dict[k]['chrom']=k
            positionlist,pedmap=variantUtils.VCF_Data.Vcf2Pedrandomdilut(**shared_dict[k],vcf_aligned=vcf_aligned)
            if options.r2:
                tn=options.outputpre+k+random_str()
                print("print in a temp map ped then invoke plink  --file  "+tn+" --allow-extra-chr --chr-set "+str(50)+" --indep-pairwise "+str(winSNPs)+" "+str(stepSNPs)+" "+options.r2+" --out "+tn)
                variantUtils.VCF_Data.printpedmap(tn,(positionlist,pedmap))
                a=0
                a+=os.system("plink --file  "+tn+" --allow-extra-chr --chr-set "+str(50)+" --indep-pairwise "+str(winSNPs)+" "+str(stepSNPs)+" "+options.r2+" --out "+tn)
                a+=os.system("plink --file  "+tn+" --allow-extra-chr --chr-set "+str(50)+" --extract "+tn+".prune.in --make-bed --recode --out "+tn+"_rfilter")
                if a!=0:
                    print("plink filter according r2 Error!")
                rest_list[allchromOrder.index(k)]=tn+"_rfilter"
                return
            try:
                rest_list[allchromOrder.index(k)]=(positionlist,pedmap)#rest_queue.put((k,(positionlist,pedmap)))#rest_list[vcfobj.chromOrder.index(k)]=
            except Exception as e:
                print(f"error catched: {e}")
                chunk_FileName=k+config.random_str()
                rest_list[allchromOrder.index(k)]=None#rest_queue.put({k:chunk_FileName})#
                pickle.dump((positionlist,pedmap),open("wapper"+chunk_FileName, 'wb'))#chunk_FileName
                rest_listError.append((k,"wapper"+chunk_FileName))
            print(shared_dict[k]["chrom"],"finished")
        def mergeMapPed(totalMapPed, mapped):
            if not totalMapPed[0] and not mapped[0]:
                return totalMapPed  # 如果都是空，直接返回
        
            totalSnpPoslist = totalMapPed[0] + mapped[0]
            for indname in mapped[1]:
                if indname in totalMapPed[1]:
                    totalMapPed[1][indname] += mapped[1][indname]
                else:
                    totalMapPed[1][indname] = mapped[1][indname]
        
            return totalSnpPoslist, totalMapPed[1]

        
        nt=int(options.nt)
        if nt >3:
            print("warnning! too many chromosomes simultaneously loaded into memory may bring some risk. make sure your machine have enough memory")
        vcfobj=variantUtils.VCF_Data(options.vcffile[0])
        vcfobj_l=[];titlelist=[];allchromOrder=[]
        for vcfFname in options.vcffile:
            vcfobj_l.append(variantUtils.VCF_Data(vcfFname))
            titlelist.append(vcfobj_l[-1].VcfIndexMap['title'])
            allchromOrder+=vcfobj_l[-1].chromOrder
        allchromOrder=list(set(allchromOrder))
        total_result = ([], {})  # 总结果是一个包含列表和字典的元组
        processes=[];maf=0
        manager = Manager()
        shared_dict = {}#manager.dict()
        rest_list=manager.list([None]*len(allchromOrder))#Queue()
        rest_listError=manager.list()
        for chromNo in allchromOrder:#vcfobj.chromOrder#for parameter in generate_parameters(vcfobj)
            #excludesitesMapBchr[chromNo] if chromNo in excludesitesMapBchr else []
            # vcfcontent=[]
            # vcfcontentl=vcfobj.getVcfListByChrom(chromNo,MQfilter=None)#vcffile.read(vcfobj.VcfIndexMap[chromNo][1]-vcffile.tell())
            excludesitesAchr=[]#chromNo = parameter["chrom"]
            # shared_dict[chromNo] = parameter    
            shared_dict[chromNo]={'vcfcontentgen':vcfobj_l,'title':titlelist,"SOFTWARE":"GATK","chrom":allchromOrder,"chromchangemap":chromchangemap,"dilutetodensity":options.dilutetodensity,"maf":float(maf),"excludesits":excludesitesAchr}
            p=Process(target=Vcf2Pedrandomdilut_Wapper,  args=(shared_dict,chromNo,rest_list,rest_listError))
            processes.append(p)
            p.start()
            if len(processes)==nt or chromNo==allchromOrder[-1]:       
                for p in processes:
                    p.join()
                processes = []
                # while not rest_queue.empty():
                #     cchrom,mapped=rest_queue.get()
                #     print("merge",cchrom)
                #     total_result=mergeMapPed(total_result, mapped)
        if rest_listError:
            for echrom,sfname in rest_listError:
                print("===========solve error chr==>>",echrom,sfname)
                l,p= pickle.load(open(sfname, 'rb'))
                total_result=mergeMapPed(total_result, (l,p))
                os.system("rm "+sfname)
        if options.r2:
            for fn in rest_list:
                if not fn or (not os.path.exists(fn+".map")):continue
                l,p=variantUtils.read_ped_map(fn)
                total_result=mergeMapPed(total_result, (l,p))
                # os.system("rm "+fn+".*")
        elif len(rest_list)-rest_list.count(None)>=1:
            total_result= reduce(mergeMapPed,[total_result]+list(rest_list))
        #positionlist,pedmap=#pool.map(Vcf2Pedrandomdilut_Wapper,parameterstuples_list)
        variantUtils.VCF_Data.printpedmap(options.outputpre,total_result)
        print("finished")
    elif  options.analysistype=="VJgendadi":
        # from multiprocessing import Pool, Lock

        outgroupidx_in_topleveltable=[8,10];minoutgroupdepth=20
        altpflagcallback=None;considerINDELandmultpleallele=False
        if not options.snpperkb  and not os.path.exists(options.outputfilename+".dilutetodensity"):
            dadisnpfN=options.outputfilename+".dilutetodensity"
        elif options.snpperkb and not os.path.exists(options.outputfilename+".dilutetodensity"+options.snpperkb.strip()):
            dadisnpfN=options.outputfilename+".dilutetodensity"+options.snpperkb.strip()
        else:
            print("file exist",os.path.exists(options.outputfilename+".dilutetodensity"))
            exit(0)
        if options.includeSNPID:
            if options.snpperkb:print("conflict -i and -d, only one allowed");exit(-1)
            includeSNPs=variantUtils.decodeSNPID(options.includeSNPID)        

        flankseqfafile=open(options.outputfilename+"."+re.search(r"[^/]*$",options.chromlist).group(0)+".fa","a")
        selectedchroms=[];chromlistfile=open(options.chromlist,"r")
        for chrrow in chromlistfile:
            chrrowlist=re.split(r'\s+',chrrow.strip())
            selectedchroms.append((chrrowlist[0].strip(),int(chrrowlist[1].strip()),int(chrrowlist[2].strip())))
        chromlistfile.close()
        import sys
        print("#",' '.join(sys.argv),file=open(options.outputfilename+'tmp','w'))

       
        totalsnp=0;totallength=0;totallengthduilt=0;totalduiltsnp=0
        """
            这里该如何多线程，同时控制好输出?一个总文件，先输出title,再来一堆个线程（chr）的输出，最终合并入总文件。
        """
        # lock = Lock()
        chridstrend=[];vcffiles=[];quantizingvcfs=[];noofindvds2quantizings=[];mode2quantizings=[];
        reffanames=[];outgroupidx_in_topleveltables=[];outpres=[];minoutgroupdepths=[];altpresflags=[];considerIndmuls=[]
        prunes=[];dilutSNPpkb=[]
        for currentchrID,currentchrstrpos,currentchrLen in selectedchroms:
            if currentchrID not in includeSNPs:continue# or int(currentchrstrpos)>includeSNPs[currentchrID][-1]
            chridstrend.append( (currentchrID,currentchrstrpos,currentchrLen))#, # fakepoolVcfmap,
            vcffiles.append(options.vcffile)
            quantizingvcfs.append(options.quantizingvcf)
            noofindvds2quantizings.append(options.noofindvds2quantizing[0])
            mode2quantizings.append(options.noofindvds2quantizing[1])
            reffanames.append(options.reffa)
            outgroupidx_in_topleveltables.append(outgroupidx_in_topleveltable) 
            outpres.append(options.outputfilename)
            minoutgroupdepths.append(minoutgroupdepth);altpresflags.append(options.alterpresence);considerIndmuls.append(considerINDELandmultpleallele)
            prunes.append((options.includeSNPID,includeSNPs))
            dilutSNPpkb.append(options.snpperkb)
        with concurrent.futures.ProcessPoolExecutor(max_workers=int(options.nt)) as executor:
            results = list(executor.map(variantUtils.GendadidsBychrom,  chridstrend,vcffiles,quantizingvcfs,noofindvds2quantizings,mode2quantizings,reffanames,outgroupidx_in_topleveltables,outpres,minoutgroupdepths,altpresflags,considerIndmuls,prunes,dilutSNPpkb) )
        # with Pool(int(options.nt)) as pool:
        #     results = pool.starmap(variantUtils.GendadidsBychrom, inputparams)
        flankseqfafile.close()#can't use in multiple processes/threads
        cmdstr="(cat "+options.outputfilename+"tmp; cat "+options.outputfilename+"."+selectedchroms[0][0]+str(selectedchroms[0][1])
        rmstr="rm "+options.outputfilename+"."+selectedchroms[0][0]+str(selectedchroms[0][1])
        for currentchrID,currentchrstrpos,currentchrLen in selectedchroms[1:]:
            #if currentchrID not in includeSNPs or int(currentchrstrpos)>includeSNPs[currentchrID][-1]:continue 
            fileiN=options.outputfilename+"."+currentchrID+str(currentchrstrpos)
            if os.path.isfile(fileiN):
                cmdstr+="; tail -n +2 "+fileiN
                rmstr=rmstr+" "+fileiN

        a=os.system(cmdstr+") > "+dadisnpfN)
        if a==0:
            print(rmstr)
            os.system(rmstr)
        for cs in results:
            print(*cs)

    elif  options.analysistype.lower()=="test":
    
        print(f"Mutect2 invoked with reference={args.reference}, input={args.input}, output={args.output}, tumor={args.tumor_sample}, normal={args.normal_sample}")
    else:
        parser.print_help()