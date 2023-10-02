
Sites = {

  'iihe' : { 
              'lsCmd'       : 'ls' ,
              'mkDir'       : False , 
              #'xrootdPath'  : 'dcap://maite.iihe.ac.be/' ,
              'xrootdPath'  : '',
              'srmPrefix'   : 'srm://maite.iihe.ac.be:8443' ,
              #'treeBaseDir' : '/pnfs/iihe/cms/store/user/xjanssen/HWWNano/' ,
              'treeBaseDir' : '/pnfs/iihe/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/' ,
              'batchQueues' : ['localgrid@cream02'],
              'slc_ver'     : 6
           } ,

  'cern' : {
              'lsCmd'       : 'ls' ,
              'mkDir'       : True ,
              #'xrootdPath'   : 'root://eoscms.cern.ch/',
              'xrootdPath'  : 'root://eosuser.cern.ch/' ,
              'xrootduserPath'  : 'root://eosuser.cern.ch/' ,
              'xrootdPath_hww'  : 'root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/',
              #'xrootduserPath_hww'  : 'root://eosuser.cern.ch//eos/user/t/tcarnaha/Summer_2022/HWW_Ntuples/',
              'treeBaseDir' : '/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/' , #I think it would be path.
              #'treeBaseDir' : '/eos/cms/store/group/phys_smp/ec/Latinos/HWWNano/',
              #'treeBaseDir':'/eos/user/t/tcarnaha/Summer_2022/HWW_Ntuples/',
               #'treeBaseDir':'/eos/user/s/sverma/Summer_2022/HWW_Ntuples/',
               'batchQueues' : ['8nh','1nd','2nd','1nw'],
               'slc_ver'     : 7
           } ,

  'sdfarm' : {
              'lsCmd'       : 'ls' ,
              'mkDir'       : False ,
	      'xrootdPath'  : 'root://cms-xrdr.private.lo:2094/',
              'treeBaseDir' : '/xrootd/store/user/jhchoi/Latino/HWWNano/',
             } ,
              #'lsCmd'       : 'xrdfs cms-xrdr.sdfarm.kr ls' ,
              #'xrootdPath'  : 'root://cms-xrdr.private.lo:2094/', inside of Korean farm
	      #'xrootdPath'  : 'root://cms-xrdr.sdfarm.kr:1094/', outside of Korean farm
              #'treeBaseDir' : '/xrd/store/user/jhchoi/Latino/HWWNano/', condor access
              #'treeBaseDir' : '/xrootd/store/user/jhchoi/Latino/HWWNano/', prompt access


  'ifca' : {
              'lsCmd'       : 'ls' ,
              'mkDir'       : True ,
              'xrootdPath'  : '' ,
              'srmPrefix'   : 'srm://srm01.ifca.es' ,
              'treeBaseDir' : '/gpfs/projects/tier3data/LatinosSkims/RunII/Nano/' ,
             },

    'kit' : {
        'lsCmd'       : 'ls' ,
        'mkDir'       : True ,
        'xrootdPath'  : '' ,
        'srmPrefix'   : 'srm://cmssrm-kit.gridka.de:8443' ,
        'treeBaseDir' : '/ceph/ntrevisa/HWWNano/' ,
        'slc_ver'     : 7
    }
}
