
import nighres
import numpy
import math
import nibabel
from PIL import Image
from nighres.io import load_volume, save_volume
import scipy.ndimage
import os
from nibabel import processing
import subprocess
import shutil
import glob

# note that all names and directories have been simplified to illustrate the process


subject='Ahead_brain_122017'

blockface = 'Ahead_brain_122017_blockface-image.nii.gz'

in_Bielschowsky = 'stains/'
in_thionin = 'stains/'
in_parvalbumin = 'stains/'
in_calbindin = 'stains/'
in_calretinin = 'stains/'

output_dir = 'coregistered/'
                
lb_Bielschowsky = 'Ahead_brain_122017_Bielschowsky_'
lb_thionin = 'Ahead_brain_122017_thionin_'
lb_parvalbumin = 'Ahead_brain_122017_parvalbumin_'
lb_calbindin = 'Ahead_brain_122017_calbindin_'
lb_calretinin = 'Ahead_brain_122017_calretinin_'

lb_blockface = subject+'bf-'

format = '.tif'

# Below are the slice numbers for each stain, along with a boolean indicating if the slice has been flipped when scanning

Bielschowsky_nums =    ['002','005','008','011','014','017','020','023','026','028','030','033','036','039','042','045','048',
                        '051','054','056','058','061','064','067','070','073','076','078','080','083','086','089','092','095','098',
                        '101','103','105','108','111','114','117','120','123','126','128','130','133','136','139','142','145','148',
                        '151','154','156','158','161','164','167','170','173','176','178','180','183','186','189','192','195','198',
                        '201','203','205','208','211','214','217','220',      '226','229',      '233',      '239',      '245',
                        '251','254',      '259',      '264',      '270','272','276','279',      '283',      '289',      '295',
                        '302',            '308',      '314',      '320',      '326','329',      '333',      '339',      '345',
                        '351',      '356','358',      '364',      '370',      '376','379',      '383',      '389',      '394',
                        '401','404',      '408',      '414',      '420',      '426','429',      '433',      '439',      '445',
                        '451','454',      '458',      '464',      '470',      '476','479',      '483','486','489','492','495','498',
                        '501','503','505','508','510','511','514','517','520','523','526','528','530','533','536','539','542','545','548',
                        '552',      '555','559','562','564','567','570','573','576','579',      '583','586','589','592','596','598',
                        '601','603','605','608','611','614','617','620','623','626','628','630','633','636','639','642','645','648',
                        '651','653','655','658','661','664','667','670','673','676','678','680','683','686','689','692','695','698',
                        '701','702','705','708','711','714','717','720','723','726','729',      '732','735','738','742','745','748',
                        '751','754','756','758','761','764','767','770','773','776','779','781','783','786','789','792','795','798',
                        '801','804','806','808','811','814','817','820','823','826','829','831','833']

Bielschowsky_orient =  [True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, 
                        True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                        True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                        True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                        True, True, True, True, True, True, True, True,       True, True,       True,       True,       True, 
                        True, True,       True,       True,       True, True, True, True,       True,       True,       True, 
                        True,             True,       True,       True,       True, True,       True,       True,       True, 
                        True,       True, True,       True,       True,       True, True,       True,       True,       True, 
                        True, True,       True,       True,       True,       True, True,       True,       True,       True, 
                        True, True,       True,       True,       True,       True, True,       True, True, True, True, True, True,
                        True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                        True,       True, True, True, True, True, True, True, True, True,       True, True, True, True, True, True,
                        True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                        True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                        True, True, True, True, True, True, True, True, True, True, True,       True, True, True, True, True, True,
                        True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                        True, True, True, True, True, True, True, True, True, True, True, True, True]

thionin_nums = ['001','004','007','010','013','016','019','022','025','029','032','035','038','041','044','047',
                '050','053','057','060','063','066','069','072','075','079','082','085','088','091','094','097',
                '100','104','107','110','113','116','119','122','125','129','132','135','138','141','144','147',
                '150','153','157','160','163','166','169','172','175','179','182','185','188','191','194','197',
                '200','204','207','210','213','216','219','222','225','228','232','235','238','241','244','247',
                '250','253','257',      '263','266','269',      '275','278','282','285','288','291','294','297',
                '300','304','307','310','313','316','319','322','325','328','332','335','338','341','344',
                '350','353','357','360','363','366','369','372','375','378','382','385','388','391','393','397', 
                '400','403','407','410','413','416','419',      '425','428','432','435','438','441','444','447',
                '450','453','457','460','463','466','469','472','475','478','482','485','488','491','494','497',
                '500','504','507',      '513','516','519','522','525','529','532','535','538','541','544','547',
                '550','554','557','561','563','566','569','572','575','578','582','585','588','591','594','597',
                '600','604','607','610','613','616','619','622','625','629','632','635','638','641','644','647',
                '650','654','657','660','663','666','669','672','675','679','682','685','688','691','694','697',
                '700','704','707','710','713','716','719','722','725','728','731','734','739','741','744','747',
                '750','753','757','760','763','766','769','772','775','778','782','785','788','791','794','797',
                '800','803','807','810','813','816','819','822','825','828','832','834']

thionin_orient=[True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True,       True, True, True,       True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, False,True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True,       True, True, True, True, True, True, True, True,
                False,True, False,True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True,       True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True]


parvalbumin_nums =      ['003','006','009','012','015','018','021','024','027','031','034','037','040','043','046','049',
                         '052','055','059','062','065','068','071','074','077','081','084','087','090','093','096','099',
                         '102','106','109','112','115','118','121','124','127','131','134','137','140','143','146','149',
                         '152','155','159','162','165','168','171','174','177','181','184','187','190','193','196','199',
                         '202','206','209','212','215','218','221','224',      '231',      '237',      '243',      '249',
                               '256',      '262',      '268',      '274',      '281',      '287',      '293',      '299',
                               '306',      '312',      '318',      '324',      '331',      '337',      '343',      '348',
                               '355',      '362',      '368',      '374',      '381',      '387',                  '399',
                               '406',      '412',      '418',      '423',      '431',      '437',      '443',      '449',
                               '456',      '462',      '468',      '474',      '481','484','487','490','493','496','499',
                         '502','506','509','512','515','518','521','524','527','531','534','537','540','543','546','549',
                         '553','556','558',      '565','568','571','574','577','581','584','587','590','593','595','599',
                         '602','606','609','612','615','618','621','624','627','631','634','637','640','643','646',
                         '652','656','659','662','665','668','671','674','677','681','684','687','690','693','696','699',
                         '703','706','709','712','715','718','721','724','727','730','733','736','740','743','746','749',
                         '752','755','759','762','765','768','771','774','777','780','784','787','790','793','796','799',
                         '802','805','809','812','815','818','821','824','827','830']
                
parvalbumin_orient =     [True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True,       True,       True,       True,       True,
                                True,       True,       True,       True,       True,       True,       True,       True,
                                True,       True,       True,       True,       True,       True,       True,       True,
                                True,       True,       True,       True,       True,       True,                   True,
                                True,       True,       True,       True,       True,       True,       True,       True,
                                True,       True,       True,       True,       True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True,       True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True,
                          True, True, True, True, True, True, True, True, True, True]

calbindin_nums =     ['227','234','240','246','252','258','265','271','277','284','290','296',
		        '303','309','315','321','327','334','340','346','352','359','365','371','377','384','390','396',
                '402','409','415','421','427','434','440','446','452','459','465','471','477']
		
calbindin_orient =   [True, True, True, True, True, True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, 
                True, True, True, True, True, True, True, True, True, True, True, True, True]

calretinin_nums  =   ['223','230','236','242','248','255','261','267','273','280','286','292','298',
                '305','311','317','323','330','336','342','347','354','361','367','373','380','386','392','398',
                '405','411','417','422','430','436','442','448','455','461','467','473','480']          
  
calretinin_orient =  [True, True, True, True, True, True, True, True, True, True, True, True, True, 
                True, True, True, True, True, True, False,False,True, True, True, True, True, True, True,
                True, True, True, True, True, True, True, True, True, True, True, True, True]
  

in_stains = [in_Bielschowsky,in_thionin,in_parvalbumin,in_calbindin,in_calretinin]
stains_nums = [Bielschowsky_nums,thionin_nums,parvalbumin_nums,calbindin_nums,calretinin_nums]          
stains_orients = [Bielschowsky_orient,thionin_orient,parvalbumin_orient,calbindin_orient,calretinin_orient]          
lb_stains = [lb_Bielschowsky,lb_thionin,lb_parvalbumin,lb_calbindin,lb_calretinin]

stains_clusters = [2,4,4,4,4]          

reverse_contrast = [True,False,True,True,True]

first = 0
last = 835
single_stain = -1

total_nums = len(calbindin_nums)+len(calretinin_nums)+len(Bielschowsky_nums)+len(thionin_nums)+len(parvalbumin_nums)

# preprocessing: scale stain and blockface images close to 200um and 
# save individual nifti images for each 2D slice

for nst,stain_images in enumerate(stains_nums):
    # offset between blockface *STACK* and staining numbers, 
    # *NOT* constant through the stack
    # (here we assume padding by 4 slices in the stack)
    bf_images = []
    bf_offset = 16
    for idx in stain_images:
        if (int(idx)<124): bf_idx = int(idx) + bf_offset
        elif (int(idx)<440): bf_idx = int(idx) + bf_offset + 1
        else: bf_idx = int(idx) + bf_offset
        bf_images.append(str(bf_idx).zfill(3))

    # resampling and cropping parameters
    stain_scale=7
    bf_scale=2 
    scaling = (0.15,0.15)

    # define here so we don't have to reload every time
    bf_img = None
        
    for idx,stain in enumerate(stain_images):
    
        # 1. Extract inverse lightness contrast & rescale it 
        output = output1_dir+lb_stains[nst]+stain+'-li-scaled-med.nii.gz'
        if (os.path.isfile(output)):
            print('1. Extract lightness: done')
            stain_nifti = output
        else:
            print('1. Extract lightness')
            # get the TIFF image
            slice_name = in_stains[nst]+lb_stains[nst]+stain+format
            if os.path.isfile(slice_name):
                slice_img = Image.open(slice_name).convert(mode='RGB')
                
                # compute simple lightness
                slice_img = numpy.array(slice_img)
                slice_li = numpy.average(slice_img,axis=2)
            
                 # reverse contrast if needed
                slice_max = numpy.max(slice_li)
                if reverse_contrast[nst]: slice_li = slice_max-slice_li
                
                # crop
                image = slice_li
                slice_li = numpy.pad(image,pad_width=((0,stain_scale),(0,stain_scale)),mode='edge')
                
                slice_crop = numpy.zeros((math.ceil(image.shape[0]/stain_scale),math.ceil(image.shape[1]/stain_scale),stain_scale*stain_scale))
                for dx in range(stain_scale):
                    for dy in range(stain_scale):
                        slice_crop[:,:,dx+dy*stain_scale] = slice_li[dx:stain_scale*math.ceil(image.shape[0]/stain_scale):stain_scale,dy:stain_scale*math.ceil(image.shape[1]/stain_scale):stain_scale]
                
                # save the median as nifti image
                sub = numpy.median(slice_crop,axis=2)
                sub = numpy.transpose(sub)
        
                header = nibabel.Nifti1Header()
                header.set_data_shape(sub.shape)
                header.set_zooms(scaling)
                
                affine = numpy.eye(4)
                affine[0,3] = -sub.shape[0]/2.0
                affine[1,3] = -sub.shape[1]/2.0
                
                # flip if needed
                if not stains_orients[nst][idx]: sub = numpy.flip(sub, axis=0)    
                
                stain_nifti = nibabel.Nifti1Image(sub,affine=affine,header=header)
                save_volume(output,stain_nifti)
                stain_nifti = output
                
            else:
                print('file not found')
                
        # 2. Extract corresponding image from pre-processed blockface
        output = output1_dir+subject+'bf-'+bf_images[idx]+'-li-inv.nii.gz'
        
        if (os.path.isfile(output)):
            print('2. Extract blockface: done')
            bf_nifti = output    
        else:
            print('2. Extract blockface')
            if (bf_img is None):
                bf_img = load_volume(blockface_img).get_fdata()
    
            bf_slice = bf_img[:,:,int(bf_images[idx])]
    
            header = nibabel.Nifti1Header()
            header.set_data_shape(bf_slice.shape)
            header.set_zooms(scaling)
    
            affine = numpy.eye(4)
            affine[0,3] = -bf_slice.shape[0]/2.0
            affine[1,3] = -bf_slice.shape[1]/2.0
            
            bf_nifti = nibabel.Nifti1Image(bf_slice,affine=affine,header=header)
            save_volume(output,bf_nifti)
        
    
# estimate the foreground and background slice by slice
print("Background estimation")

if (first==0 and len(glob.glob(output2_dir+'foreground-[0-9][0-9][0-9].nii.gz'))==total_nums):
    print(" fully completed")
else:
    for num in range(first,last+1):
        num = str(int(num)).zfill(3)
        for idx,stain in enumerate(stains_nums):
            if (num in stain and (single_stain<0 or idx==single_stain) ):
                output = output2_dir+'foreground-'+num+'.nii.gz'
                if (os.path.isfile(output)):
                    print('Background estimation: done')
                else:
                    histo = output1_dir+lb_stains[idx]+num+'-li-scaled-med.nii.gz'
                    histon4 = output2_dir+lb_stains[idx]+num+'-li-scaled-med-n4.nii.gz'
    
                    # N4 inhomogeneity correction
                    command = 'N4BiasFieldCorrection -d 2 -i '+histo+' -o '+histon4
                    print(command)
                    try:
                        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
                    except subprocess.CalledProcessError as e:
                        msg = 'execution failed (error code '+str(e.returncode)+')\n Output: '+str(e.output)
                        raise subprocess.CalledProcessError(msg)
    
                    # Fuzzy C-means segmentation
                    fcm = nighres.segmentation.fuzzy_cmeans(histon4, clusters=stains_clusters[idx], max_iterations=50, max_difference=0.01, 
                        smoothing=0.0, fuzziness=2.0, mask_zero=True,
                        save_data=False)
                    
                    mask_data = fcm['classification'].get_fdata()==0
                    bg_data = fcm['memberships'][0].get_fdata()+mask_data
                    
                    bg = nibabel.Nifti1Image(bg_data, None, None)
    
                    # intensity propagation
                    propag = nighres.intensity.intensity_propagation(bg, mask=None, combine='mean', distance_mm=10.0,
                          target='lower', scaling=0.9, save_data=False)
                    
                    fg = nibabel.Nifti1Image(1.0-propag['result'].get_fdata(), None, None)
                    nighres.io.save_volume(output, fg)



print("First multi-stain registration: forward")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov01-*_ants-def0.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for num in range(first,last+1):
        num = str(int(num)).zfill(3)
        for idx,stain in enumerate(stains_nums):
            if num in stain:
                #histo = output1_dir+lb_stains[idx]+num+'-li-scaled-med.nii.gz'
                histo = output2_dir+lb_stains[idx]+num+'-li-scaled-med-n4.nii.gz'
                histo_vessel = output2_dir+'foreground-'+num+'.nii.gz'
                print("current: "+histo_vessel)
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)

                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                
                # previous slices
                histo_prev10 = None
                histo_prev11 = None
                # try num-1
                prev1 = str(int(num)-1).zfill(3)
                for pid,prev_stain in enumerate(stains_nums):
                    if prev1 in prev_stain:
                        histo_prev10 = output3_dir+'foreground-markov01-'+prev1+'_ants-def0.nii.gz'
                        histo_prev11 = output3_dir+'foreground-markov01-'+prev1+'_ants-def1.nii.gz'
                        print("first previous: "+histo_prev10)
                        print("first previous: "+histo_prev11)
        
                histo_prev20 = None
                histo_prev21 = None
                # try num-2
                prev2 = str(int(num)-2).zfill(3)
                for pid,prev_stain in enumerate(stains_nums):
                    if prev2 in prev_stain:
                        histo_prev20 = output3_dir+'foreground-markov01-'+prev2+'_ants-def0.nii.gz'
                        histo_prev21 = output3_dir+'foreground-markov01-'+prev2+'_ants-def1.nii.gz'
                        print("second previous: "+histo_prev20)
                        print("second previous: "+histo_prev21)
        
                # look for previous from same stain, use original contrast
                histo_prev3 = None
                loc = stain.index(num)
                if loc>0:
                    histo_prev3 = output3_dir+'foreground-markov01-'+stain[loc-1]+'_ants-def0.nii.gz'
                    print("previous (same contrast): "+histo_prev3)
                    
                sources = [histo]
                targets = [blockface]
                image_weights = [1.0]
                
                # start with a vessel map in place 1
                if histo_prev11 is not None and os.path.exists(histo_prev11):
                    sources.append(histo_vessel)
                    targets.append(histo_prev11)
                    image_weights.append(0.25)
                     
                if histo_prev10 is not None and os.path.exists(histo_prev10):
                    sources.append(histo)
                    targets.append(histo_prev10)
                    image_weights.append(0.25)
                     
                if histo_prev21 is not None and os.path.exists(histo_prev21):
                    sources.append(histo_vessel)
                    targets.append(histo_prev21)
                    image_weights.append(0.25)
                    
                if histo_prev20 is not None and os.path.exists(histo_prev20):
                    sources.append(histo)
                    targets.append(histo_prev20)
                    image_weights.append(0.25)
                    
                if histo_prev3 is not None and os.path.exists(histo_prev3):
                    sources.append(histo)
                    targets.append(histo_prev3)
                    image_weights.append(0.5)
                    
                if len(sources)==1:
                    sources.append(histo_vessel)
                    targets.append(blockface)
                    image_weights.append(0.5)
                                        
                nighres.registration.embedded_antsreg_2d_multi(source_images=sources, 
                            target_images=targets, image_weights=image_weights,
                            run_rigid=True,
                            rigid_iterations=1000,
                            run_affine=True,
                            affine_iterations=2000,
                            run_syn=True,
                            coarse_iterations=2000,
                            medium_iterations=1000, fine_iterations=200,
                            scaling_factor=64,
                            cost_function='MutualInformation',
                            interpolation='NearestNeighbor',
                            regularization='High',
                            convergence=1e-6,
                            mask_zero=False,
                            ignore_affine=True, ignore_orient=True, ignore_res=True,
                            save_data=True, overwrite=(first!=0),
                            output_dir=output3_dir,
                            file_name='foreground-markov01-'+num)


# second pass based on co-aligned results: backwards
print("First co-registration: backward")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov02-*_ants-def0.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for num in range(last,first-1,-1):
        num = str(int(num)).zfill(3)
        for idx,stain in enumerate(stains_nums):
            if num in stain:
                histo = output2_dir+lb_stains[idx]+num+'-li-scaled-med-n4.nii.gz'
                histo_vessel = output2_dir+'foreground-'+num+'.nii.gz'
                print("current: "+histo_vessel)
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                print("ref: "+blockface)
                
                # next slices
                histo_next10 = None
                histo_next11 = None
                # try num+1
                next1 = str(int(num)+1).zfill(3)
                for nid,next_stain in enumerate(stains_nums):
                    if next1 in next_stain:
                        histo_next10 = output3_dir+'foreground-markov02-'+next1+'_ants-def0.nii.gz'
                        histo_next11 = output3_dir+'foreground-markov02-'+next1+'_ants-def1.nii.gz'
                        print("first next: "+histo_next10)
                        print("first next: "+histo_next11)
                
                histo_next20 = None
                histo_next21 = None
                # try num+2
                next2 = str(int(num)+2).zfill(3)
                for nid,next_stain in enumerate(stains_nums):
                    if next2 in next_stain:
                        histo_next20 = output3_dir+'foreground-markov02-'+next2+'_ants-def0.nii.gz'
                        histo_next21 = output3_dir+'foreground-markov02-'+next2+'_ants-def1.nii.gz'
                        print("second next: "+histo_next20)
                        print("second next: "+histo_next21)
    
                # look for previous from same stain, use original contrast
                histo_next3 = None
                loc = stain.index(num)
                if loc<len(stain)-1:
                    histo_next3 = output3_dir+'foreground-markov02-'+stain[loc+1]+'_ants-def0.nii.gz'
                    print("previous (same contrast): "+histo_next3)
                    
                sources = [histo]
                targets = [blockface]
                image_weights = [1.0]
                
                # start with a vessel map in place 1
                if histo_next11 is not None and os.path.exists(histo_next11):
                    sources.append(histo_vessel)
                    targets.append(histo_next11)
                    image_weights.append(0.25)
                    
                if histo_next10 is not None and os.path.exists(histo_next10):
                    sources.append(histo)
                    targets.append(histo_next10)
                    image_weights.append(0.25)
                    
                if histo_next21 is not None and os.path.exists(histo_next21):
                    sources.append(histo_vessel)
                    targets.append(histo_next21)
                    image_weights.append(0.25)
                    
                if histo_next20 is not None and os.path.exists(histo_next20):
                    sources.append(histo)
                    targets.append(histo_next20)
                    image_weights.append(0.25)
                    
                if histo_next3 is not None and os.path.exists(histo_next3):
                    sources.append(histo)
                    targets.append(histo_next3)
                    image_weights.append(0.5)
                    
                if len(sources)==1:
                    sources.append(histo_vessel)
                    targets.append(blockface)
                    image_weights.append(0.5)
                    
                nighres.registration.embedded_antsreg_2d_multi(source_images=sources, 
                            target_images=targets, image_weights=image_weights,
                            run_rigid=True,
                            rigid_iterations=1000,
                            run_affine=True,
                            affine_iterations=2000,
                            run_syn=True,
                            coarse_iterations=2000,
                            medium_iterations=1000, fine_iterations=200,
                            scaling_factor=64,
                            cost_function='MutualInformation',
                            interpolation='NearestNeighbor',
                            regularization='High',
                            convergence=1e-6,
                            mask_zero=False,
                            ignore_affine=True, ignore_orient=True, ignore_res=True,
                            save_data=True, overwrite=(first!=0),
                            output_dir=output3_dir,
                            file_name='foreground-markov02-'+num)

# keep best result
print("Select best result")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov12-*_ants-invmap.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for idx,stain in enumerate(stains_nums):
        for num in stain:
            if int(num)>=first and int(num)<=last:
            
                # select the highest MI result for output
                histo1 = output3_dir+'foreground-markov01-'+num+'_ants-def0.nii.gz'
                vessel1 = output3_dir+'foreground-markov01-'+num+'_ants-def1.nii.gz'
                mapping1 = output3_dir+'foreground-markov01-'+num+'_ants-map.nii.gz'
                inverse1 = output3_dir+'foreground-markov01-'+num+'_ants-invmap.nii.gz'
                
                histo2 = output3_dir+'foreground-markov02-'+num+'_ants-def0.nii.gz'
                vessel2 = output3_dir+'foreground-markov02-'+num+'_ants-def1.nii.gz'
                mapping2 = output3_dir+'foreground-markov02-'+num+'_ants-map.nii.gz'
                inverse2 = output3_dir+'foreground-markov02-'+num+'_ants-invmap.nii.gz'
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                
                curr1 = nighres.io.load_volume(histo1).get_fdata()
                curr2 = nighres.io.load_volume(histo2).get_fdata()
                curr = nighres.io.load_volume(blockface).get_fdata()
                    
                nonzero = curr>0
                curr1 = curr1[nonzero]
                curr2 = curr2[nonzero]
                curr = curr[nonzero]  
    
                p1,v1 = numpy.histogram(curr1.flatten(), bins=100, density=True)
                p2,v2 = numpy.histogram(curr2.flatten(), bins=100, density=True)
                pc,vc = numpy.histogram(curr.flatten(), bins=100, density=True)
                    
                p1c,v1,vc = numpy.histogram2d(curr1.flatten(), curr.flatten(), bins=100, density=True)
                p2c,v2,vc = numpy.histogram2d(curr2.flatten(), curr.flatten(), bins=100, density=True)
                
                p1pc = numpy.outer(p1,pc)
                p2pc = numpy.outer(p2,pc)
                         
                mi1c = numpy.sum(p1c*numpy.log(p1c/(p1pc),where=(p1c*p1pc>0)))
                mi2c = numpy.sum(p2c*numpy.log(p2c/(p2pc),where=(p2c*p2pc>0)))
                
                print("MI: "+str(mi1c)+", "+str(mi2c)+" -> "+str(mi2c>mi1c))
                
                target = output3_dir+'foreground-markov12-'+num+'_ants-def0.nii.gz'
                vessel = output3_dir+'foreground-markov12-'+num+'_ants-def1.nii.gz'
                mapping = output3_dir+'foreground-markov12-'+num+'_ants-map.nii.gz'
                inverse = output3_dir+'foreground-markov12-'+num+'_ants-invmap.nii.gz'
                if (mi1c>mi2c):
                    shutil.copyfile(histo1,target)
                    shutil.copyfile(vessel1,vessel)
                    shutil.copyfile(mapping1,mapping)
                    shutil.copyfile(inverse1,inverse)
                else:
                    shutil.copyfile(histo2,target)
                    shutil.copyfile(vessel2,vessel)
                    shutil.copyfile(mapping2,mapping)
                    shutil.copyfile(inverse2,inverse)
        

print("Second multi-stain registration: forward")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov3-*_ants-def0.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for num in range(first,last+1):
        num = str(int(num)).zfill(3)
        for idx,stain in enumerate(stains_nums):
            if num in stain:
                histo = output3_dir+'foreground-markov12-'+num+'_ants-def0.nii.gz'
                histo_vessel = output3_dir+'foreground-markov12-'+num+'_ants-def1.nii.gz'
                print("current: "+histo)
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                print("ref: "+blockface)
                
                # previous slices
                histo_prev10 = None
                histo_prev11 = None
                # try num-1
                prev1 = str(int(num)-1).zfill(3)
                for pid,prev_stain in enumerate(stains_nums):
                    if prev1 in prev_stain:
                        histo_prev10 = output3_dir+'foreground-markov3-'+prev1+'_ants-def0.nii.gz'
                        histo_prev11 = output3_dir+'foreground-markov3-'+prev1+'_ants-def1.nii.gz'
                        print("first previous: "+histo_prev10)
                        print("first previous: "+histo_prev11)
        
                histo_prev20 = None
                histo_prev21 = None
                # try num-2
                prev2 = str(int(num)-2).zfill(3)
                for pid,prev_stain in enumerate(stains_nums):
                    if prev2 in prev_stain:
                        histo_prev20 = output3_dir+'foreground-markov3-'+prev2+'_ants-def0.nii.gz'
                        histo_prev21 = output3_dir+'foreground-markov3-'+prev2+'_ants-def1.nii.gz'
                        print("second previous: "+histo_prev20)
                        print("second previous: "+histo_prev21)
        
                # look for previous from same stain, use original contrast
                histo_prev3 = None
                loc = stain.index(num)
                if loc>0:
                    histo_prev3 = output3_dir+'foreground-markov3-'+stain[loc-1]+'_ants-def0.nii.gz'
                    print("previous (same contrast): "+histo_prev3)
                    
                sources = [histo]
                targets = [blockface]
                image_weights = [1.0]
                
                # start with a vessel map in place 1
                if histo_prev11 is not None and os.path.exists(histo_prev11):
                    sources.append(histo_vessel)
                    targets.append(histo_prev11)
                    image_weights.append(0.25)
                     
                if histo_prev10 is not None and os.path.exists(histo_prev10):
                    sources.append(histo)
                    targets.append(histo_prev10)
                    image_weights.append(0.25)
                     
                if histo_prev21 is not None and os.path.exists(histo_prev21):
                    sources.append(histo_vessel)
                    targets.append(histo_prev21)
                    image_weights.append(0.25)
                    
                if histo_prev20 is not None and os.path.exists(histo_prev20):
                    sources.append(histo)
                    targets.append(histo_prev20)
                    image_weights.append(0.25)
                    
                if histo_prev3 is not None and os.path.exists(histo_prev3):
                    sources.append(histo)
                    targets.append(histo_prev3)
                    image_weights.append(0.5)
                    
                if len(sources)==1:
                    sources.append(histo_vessel)
                    targets.append(blockface)
                    image_weights.append(0.5)
                                        
                nighres.registration.embedded_antsreg_2d_multi(source_images=sources, 
                            target_images=targets, image_weights=image_weights,
                            run_rigid=True,
                            rigid_iterations=1000,
                            run_affine=True,
                            affine_iterations=2000,
                            run_syn=True,
                            coarse_iterations=2000,
                            medium_iterations=1000, fine_iterations=200,
                            scaling_factor=64,
                            cost_function='MutualInformation',
                            interpolation='NearestNeighbor',
                            regularization='High',
                            convergence=1e-6,
                            mask_zero=False,
                            ignore_affine=True, ignore_orient=True, ignore_res=True,
                            save_data=True, overwrite=(first!=0),
                            output_dir=output3_dir,
                            file_name='foreground-markov3-'+num)



# keep best result
print("Select best result")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov123-*_ants-invmap.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for idx,stain in enumerate(stains_nums):
        for num in stain:
            if int(num)>=first and int(num)<=last:
            
                # select the highest MI result for output
                histo12 = output3_dir+'foreground-markov12-'+num+'_ants-def0.nii.gz'
                vessel12 = output3_dir+'foreground-markov12-'+num+'_ants-def1.nii.gz'
                mapping12 = output3_dir+'foreground-markov12-'+num+'_ants-map.nii.gz'
                inverse12 = output3_dir+'foreground-markov12-'+num+'_ants-invmap.nii.gz'
                
                histo3 = output3_dir+'foreground-markov3-'+num+'_ants-def0.nii.gz'
                vessel3 = output3_dir+'foreground-markov3-'+num+'_ants-def1.nii.gz'
                mapping3 = output3_dir+'foreground-markov3-'+num+'_ants-map.nii.gz'
                inverse3 = output3_dir+'foreground-markov3-'+num+'_ants-invmap.nii.gz'
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                
                curr12 = nighres.io.load_volume(histo12).get_fdata()
                curr3 = nighres.io.load_volume(histo3).get_fdata()
                curr = nighres.io.load_volume(blockface).get_fdata()
                    
                nonzero = curr>0
                curr12 = curr12[nonzero]
                curr3 = curr3[nonzero]
                curr = curr[nonzero]  
    
                p12,v12 = numpy.histogram(curr12.flatten(), bins=100, density=True)
                p3,v3 = numpy.histogram(curr3.flatten(), bins=100, density=True)
                pc,vc = numpy.histogram(curr.flatten(), bins=100, density=True)
                    
                p12c,v12,vc = numpy.histogram2d(curr12.flatten(), curr.flatten(), bins=100, density=True)
                p3c,v3,vc = numpy.histogram2d(curr3.flatten(), curr.flatten(), bins=100, density=True)
                
                p12pc = numpy.outer(p12,pc)
                p3pc = numpy.outer(p3,pc)
                         
                mi12c = numpy.sum(p12c*numpy.log(p12c/(p12pc),where=(p12c*p12pc>0)))
                mi3c = numpy.sum(p3c*numpy.log(p3c/(p3pc),where=(p3c*p3pc>0)))
                
                print("MI: "+str(mi12c)+", "+str(mi3c)+" -> "+str(mi3c>mi12c))
                
                target = output3_dir+'foreground-markov123-'+num+'_ants-def0.nii.gz'
                vessel = output3_dir+'foreground-markov123-'+num+'_ants-def1.nii.gz'
                mapping = output3_dir+'foreground-markov123-'+num+'_ants-map.nii.gz'
                inverse = output3_dir+'foreground-markov123-'+num+'_ants-invmap.nii.gz'
                if (mi12c>mi3c):
                    shutil.copyfile(histo12,target)
                    shutil.copyfile(vessel12,vessel)
                    shutil.copyfile(mapping12,mapping)
                    shutil.copyfile(inverse12,inverse)
                else:
                    shutil.copyfile(histo3,target)
                    shutil.copyfile(vessel3,vessel)
                    combined = nighres.registration.apply_coordinate_mappings_2d(image=mapping12, mapping1=mapping3,
                                interpolation="nearest", padding="zero",
                                save_data=False)

                    nighres.io.save_volume(mapping,combined['result'])
                
                    combined = nighres.registration.apply_coordinate_mappings_2d(image=inverse3, mapping1=inverse12,
                                interpolation="nearest", padding="zero",
                                save_data=False)

                    nighres.io.save_volume(inverse,combined['result'])
                

# second pass based on co-aligned results: backwards
print("First co-registration: backward")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov4-*_ants-def0.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for num in range(last,first-1,-1):
        num = str(int(num)).zfill(3)
        for idx,stain in enumerate(stains_nums):
            if num in stain:
                histo = output3_dir+'foreground-markov123-'+num+'_ants-def0.nii.gz'
                histo_vessel = output3_dir+'foreground-markov123-'+num+'_ants-def1.nii.gz'
                print("current: "+histo)
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                print("ref: "+blockface)
                
                # next slices
                histo_next10 = None
                histo_next11 = None
                # try num+1
                next1 = str(int(num)+1).zfill(3)
                for nid,next_stain in enumerate(stains_nums):
                    if next1 in next_stain:
                        histo_next10 = output3_dir+'foreground-markov4-'+next1+'_ants-def0.nii.gz'
                        histo_next11 = output3_dir+'foreground-markov4-'+next1+'_ants-def1.nii.gz'
                        print("first next: "+histo_next10)
                        print("first next: "+histo_next11)
                
                histo_next20 = None
                histo_next21 = None
                # try num+2
                next2 = str(int(num)+2).zfill(3)
                for nid,next_stain in enumerate(stains_nums):
                    if next2 in next_stain:
                        histo_next20 = output3_dir+'foreground-markov4-'+next2+'_ants-def0.nii.gz'
                        histo_next21 = output3_dir+'foreground-markov4-'+next2+'_ants-def1.nii.gz'
                        print("second next: "+histo_next20)
                        print("second next: "+histo_next21)
    
                # look for previous from same stain, use original contrast
                histo_next3 = None
                loc = stain.index(num)
                if loc<len(stain)-1:
                    histo_next3 = output3_dir+'foreground-markov4-'+stain[loc+1]+'_ants-def0.nii.gz'
                    print("previous (same contrast): "+histo_next3)
                    
                sources = [histo]
                targets = [blockface]
                image_weights = [1.0]
                
                # start with a vessel map in place 1
                if histo_next11 is not None and os.path.exists(histo_next11):
                    sources.append(histo_vessel)
                    targets.append(histo_next11)
                    image_weights.append(0.25)
                    
                if histo_next10 is not None and os.path.exists(histo_next10):
                    sources.append(histo)
                    targets.append(histo_next10)
                    image_weights.append(0.25)
                    
                if histo_next21 is not None and os.path.exists(histo_next21):
                    sources.append(histo_vessel)
                    targets.append(histo_next21)
                    image_weights.append(0.25)
                    
                if histo_next20 is not None and os.path.exists(histo_next20):
                    sources.append(histo)
                    targets.append(histo_next20)
                    image_weights.append(0.25)
                    
                if histo_next3 is not None and os.path.exists(histo_next3):
                    sources.append(histo)
                    targets.append(histo_next3)
                    image_weights.append(0.5)
                    
                if len(sources)==1:
                    sources.append(histo_vessel)
                    targets.append(blockface)
                    image_weights.append(0.5)
                    
                nighres.registration.embedded_antsreg_2d_multi(source_images=sources, 
                            target_images=targets, image_weights=image_weights,
                            run_rigid=False,
                            rigid_iterations=1000,
                            run_affine=False,
                            affine_iterations=2000,
                            run_syn=True,
                            coarse_iterations=2000,
                            medium_iterations=1000, fine_iterations=200,
                            scaling_factor=64,
                            cost_function='MutualInformation',
                            interpolation='NearestNeighbor',
                            regularization='High',
                            convergence=1e-6,
                            mask_zero=False,
                            ignore_affine=True, ignore_orient=True, ignore_res=True,
                            save_data=True, overwrite=(first!=0),
                            output_dir=output3_dir,
                            file_name='foreground-markov4-'+num)


# keep best result
print("Select best result")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov1234-*_ants-invmap.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for idx,stain in enumerate(stains_nums):
        for num in stain:
            if int(num)>=first and int(num)<=last:
            
                # select the highest MI result for output
                histo123 = output3_dir+'foreground-markov123-'+num+'_ants-def0.nii.gz'
                vessel123 = output3_dir+'foreground-markov123-'+num+'_ants-def1.nii.gz'
                mapping123 = output3_dir+'foreground-markov123-'+num+'_ants-map.nii.gz'
                inverse123 = output3_dir+'foreground-markov123-'+num+'_ants-invmap.nii.gz'
                
                histo4 = output3_dir+'foreground-markov4-'+num+'_ants-def0.nii.gz'
                vessel4 = output3_dir+'foreground-markov4-'+num+'_ants-def1.nii.gz'
                mapping4 = output3_dir+'foreground-markov4-'+num+'_ants-map.nii.gz'
                inverse4 = output3_dir+'foreground-markov4-'+num+'_ants-invmap.nii.gz'
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                
                curr123 = nighres.io.load_volume(histo123).get_fdata()
                curr4 = nighres.io.load_volume(histo4).get_fdata()
                curr = nighres.io.load_volume(blockface).get_fdata()
                    
                nonzero = curr>0
                curr123 = curr123[nonzero]
                curr4 = curr4[nonzero]
                curr = curr[nonzero]  
    
                p123,v123 = numpy.histogram(curr123.flatten(), bins=100, density=True)
                p4,v4 = numpy.histogram(curr4.flatten(), bins=100, density=True)
                pc,vc = numpy.histogram(curr.flatten(), bins=100, density=True)
                    
                p123c,v123,vc = numpy.histogram2d(curr123.flatten(), curr.flatten(), bins=100, density=True)
                p4c,v4,vc = numpy.histogram2d(curr4.flatten(), curr.flatten(), bins=100, density=True)
                
                p123pc = numpy.outer(p123,pc)
                p4pc = numpy.outer(p4,pc)
                         
                mi123c = numpy.sum(p123c*numpy.log(p123c/(p123pc),where=(p123c*p123pc>0)))
                mi4c = numpy.sum(p4c*numpy.log(p4c/(p4pc),where=(p4c*p4pc>0)))
                
                print("MI: "+str(mi123c)+", "+str(mi4c)+" -> "+str(mi4c>mi123c))
                
                target = output3_dir+'foreground-markov1234-'+num+'_ants-def0.nii.gz'
                vessel = output3_dir+'foreground-markov1234-'+num+'_ants-def1.nii.gz'
                mapping = output3_dir+'foreground-markov1234-'+num+'_ants-map.nii.gz'
                inverse = output3_dir+'foreground-markov1234-'+num+'_ants-invmap.nii.gz'
                if (mi123c>mi4c):
                    shutil.copyfile(histo123,target)
                    shutil.copyfile(vessel123,vessel)
                    shutil.copyfile(mapping123,mapping)
                    shutil.copyfile(inverse123,inverse)
                elif (mi4c>mi123c):
                    shutil.copyfile(histo4,target)
                    shutil.copyfile(vessel4,vessel)
                    combined = nighres.registration.apply_coordinate_mappings_2d(image=mapping123, mapping1=mapping4,
                                interpolation="nearest", padding="zero",
                                save_data=False)

                    nighres.io.save_volume(mapping,combined['result'])
                
                    combined = nighres.registration.apply_coordinate_mappings_2d(image=inverse4, mapping1=inverse123,
                                interpolation="nearest", padding="zero",
                                save_data=False)

                    nighres.io.save_volume(inverse,combined['result'])
                
        
# third pass based on co-aligned results one last time: forward
print("Third co-registration: forward")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov5-*_ants-def0.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for num in range(first,last+1):
        num = str(int(num)).zfill(3)
        for idx,stain in enumerate(stains_nums):
            if num in stain:
                histo = output3_dir+'foreground-markov1234-'+num+'_ants-def0.nii.gz'
                histo_vessel = output3_dir+'foreground-markov1234-'+num+'_ants-def1.nii.gz'
                print("current: "+histo)
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                print("ref: "+blockface)
                
                # previous slices
                histo_prev10 = None
                histo_prev11 = None
                # try num-1
                prev1 = str(int(num)-1).zfill(3)
                for pid,prev_stain in enumerate(stains_nums):
                    if prev1 in prev_stain:
                        histo_prev10 = output3_dir+'foreground-markov5-'+prev1+'_ants-def0.nii.gz'
                        histo_prev11 = output3_dir+'foreground-markov5-'+prev1+'_ants-def1.nii.gz'
                        print("first previous: "+histo_prev10)
                        print("first previous: "+histo_prev11)
        
                histo_prev20 = None
                histo_prev21 = None
                # try num-2
                prev2 = str(int(num)-2).zfill(3)
                for pid,prev_stain in enumerate(stains_nums):
                    if prev2 in prev_stain:
                        histo_prev20 = output3_dir+'foreground-markov5-'+prev2+'_ants-def0.nii.gz'
                        histo_prev21 = output3_dir+'foreground-markov5-'+prev2+'_ants-def1.nii.gz'
                        print("second previous: "+histo_prev20)
                        print("second previous: "+histo_prev21)
        
                # look for previous from same stain, use original contrast
                histo_prev3 = None
                loc = stain.index(num)
                if loc>0:
                    histo_prev3 = output3_dir+'foreground-markov5-'+stain[loc-1]+'_ants-def0.nii.gz'
                    print("previous (same contrast): "+histo_prev3)
                    
                sources = [histo]
                targets = [blockface]
                image_weights = [1.0]
                
                # start with a vessel map in place 1
                if histo_prev11 is not None and os.path.exists(histo_prev11):
                    sources.append(histo_vessel)
                    targets.append(histo_prev11)
                    image_weights.append(0.25)
                     
                if histo_prev10 is not None and os.path.exists(histo_prev10):
                    sources.append(histo)
                    targets.append(histo_prev10)
                    image_weights.append(0.25)
                     
                if histo_prev21 is not None and os.path.exists(histo_prev21):
                    sources.append(histo_vessel)
                    targets.append(histo_prev21)
                    image_weights.append(0.25)
                    
                if histo_prev20 is not None and os.path.exists(histo_prev20):
                    sources.append(histo)
                    targets.append(histo_prev20)
                    image_weights.append(0.25)
                    
                if histo_prev3 is not None and os.path.exists(histo_prev3):
                    sources.append(histo)
                    targets.append(histo_prev3)
                    image_weights.append(0.5)
                    
                if len(sources)==1:
                    sources.append(histo_vessel)
                    targets.append(blockface)
                    image_weights.append(0.5)
                                        
                nighres.registration.embedded_antsreg_2d_multi(source_images=sources, 
                            target_images=targets, image_weights=image_weights,
                            run_rigid=False,
                            rigid_iterations=1000,
                            run_affine=False,
                            affine_iterations=2000,
                            run_syn=True,
                            coarse_iterations=2000,
                            medium_iterations=1000, fine_iterations=200,
                            scaling_factor=64,
                            cost_function='MutualInformation',
                            interpolation='NearestNeighbor',
                            regularization='High',
                            convergence=1e-6,
                            mask_zero=False,
                            ignore_affine=True, ignore_orient=True, ignore_res=True,
                            save_data=True, overwrite=(first!=0),
                            output_dir=output3_dir,
                            file_name='foreground-markov5-'+num)

print("Third co-registration: backward")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov6-*_ants-def0.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for num in range(last,first-1,-1):
        num = str(int(num)).zfill(3)
        for idx,stain in enumerate(stains_nums):
            if num in stain:
                histo = output3_dir+'foreground-markov5-'+num+'_ants-def0.nii.gz'
                histo_vessel = output3_dir+'foreground-markov5-'+num+'_ants-def1.nii.gz'
                print("current: "+histo)
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                print("ref: "+blockface)
                
                # next slices
                histo_next10 = None
                histo_next11 = None
                # try num+1
                next1 = str(int(num)+1).zfill(3)
                for nid,next_stain in enumerate(stains_nums):
                    if next1 in next_stain:
                        histo_next10 = output3_dir+'foreground-markov6-'+next1+'_ants-def0.nii.gz'
                        histo_next11 = output3_dir+'foreground-markov6-'+next1+'_ants-def1.nii.gz'
                        print("first next: "+histo_next10)
                        print("first next: "+histo_next11)
                
                histo_next20 = None
                histo_next21 = None
                # try num+2
                next2 = str(int(num)+2).zfill(3)
                for nid,next_stain in enumerate(stains_nums):
                    if next2 in next_stain:
                        histo_next20 = output3_dir+'foreground-markov6-'+next2+'_ants-def0.nii.gz'
                        histo_next21 = output3_dir+'foreground-markov6-'+next2+'_ants-def1.nii.gz'
                        print("second next: "+histo_next20)
                        print("second next: "+histo_next21)
    
                # look for previous from same stain, use original contrast
                histo_next3 = None
                loc = stain.index(num)
                if loc<len(stain)-1:
                    histo_next3 = output3_dir+'foreground-markov6-'+stain[loc+1]+'_ants-def0.nii.gz'
                    print("previous (same contrast): "+histo_next3)
                    
                sources = [histo]
                targets = [blockface]
                image_weights = [1.0]
                
                # start with a vessel map in place 1
                if histo_next11 is not None and os.path.exists(histo_next11):
                    sources.append(histo_vessel)
                    targets.append(histo_next11)
                    image_weights.append(0.25)
                    
                if histo_next10 is not None and os.path.exists(histo_next10):
                    sources.append(histo)
                    targets.append(histo_next10)
                    image_weights.append(0.25)
                    
                if histo_next21 is not None and os.path.exists(histo_next21):
                    sources.append(histo_vessel)
                    targets.append(histo_next21)
                    image_weights.append(0.25)
                    
                if histo_next20 is not None and os.path.exists(histo_next20):
                    sources.append(histo)
                    targets.append(histo_next20)
                    image_weights.append(0.25)
                    
                if histo_next3 is not None and os.path.exists(histo_next3):
                    sources.append(histo)
                    targets.append(histo_next3)
                    image_weights.append(0.5)
                    
                if len(sources)==1:
                    sources.append(histo_vessel)
                    targets.append(blockface)
                    image_weights.append(0.5)
                    
                nighres.registration.embedded_antsreg_2d_multi(source_images=sources, 
                            target_images=targets, image_weights=image_weights,
                            run_rigid=False,
                            rigid_iterations=1000,
                            run_affine=False,
                            affine_iterations=2000,
                            run_syn=True,
                            coarse_iterations=2000,
                            medium_iterations=1000, fine_iterations=200,
                            scaling_factor=64,
                            cost_function='MutualInformation',
                            interpolation='NearestNeighbor',
                            regularization='High',
                            convergence=1e-6,
                            mask_zero=False,
                            ignore_affine=True, ignore_orient=True, ignore_res=True,
                            save_data=True, overwrite=(first!=0),
                            output_dir=output3_dir,
                            file_name='foreground-markov6-'+num)


# combine mappings to interpolate original stains
print("Use final result (for best smoothness)")
if (first==0 and len(glob.glob(output3_dir+'foreground-markov-final-*_ants-invmap.nii.gz'))==total_nums):
    print(" fully completed")
else:
    for idx,stain in enumerate(stains_nums):
        for num in stain:
            if int(num)>=first and int(num)<=last:
            
                # select the highest MI result for output: not here
                mapping1234 = output3_dir+'foreground-markov1234-'+num+'_ants-map.nii.gz'
                inverse1234 = output3_dir+'foreground-markov1234-'+num+'_ants-invmap.nii.gz'
                
                mapping5 = output3_dir+'foreground-markov5-'+num+'_ants-map.nii.gz'
                inverse5 = output3_dir+'foreground-markov5-'+num+'_ants-invmap.nii.gz'
                
                mapping6 = output3_dir+'foreground-markov6-'+num+'_ants-map.nii.gz'
                inverse6 = output3_dir+'foreground-markov6-'+num+'_ants-invmap.nii.gz'
                
                mapping = output3_dir+'foreground-markov-final-'+num+'_ants-map.nii.gz'
                inverse = output3_dir+'foreground-markov-final-'+num+'_ants-invmap.nii.gz'
 
                combined = nighres.registration.apply_coordinate_mappings_2d(image=mapping1234, 
                                mapping1=mapping5, mapping2=mapping6,
                                interpolation="nearest", padding="zero",
                                save_data=False)

                nighres.io.save_volume(mapping,combined['result'])
                
                combined = nighres.registration.apply_coordinate_mappings_2d(image=inverse6, 
                                mapping1=inverse5, mapping2=inverse1234,
                                interpolation="nearest", padding="zero",
                                save_data=False)

                nighres.io.save_volume(inverse,combined['result'])
                
                histo = output2_dir+lb_stains[idx]+num+'-li-scaled-med-n4.nii.gz'
                nighres.registration.apply_coordinate_mappings_2d(image=histo, mapping1=mapping,
                                interpolation="nearest", padding="zero",
                                save_data=True, overwrite=(first!=0), 
                                output_dir=output3_dir,
                                file_name='foreground-markov-final-'+num+'.nii.gz')



# stack all foreground in bf space
print("Stack best result")
if (first==0 and os.path.isfile(output3_dir+"foreground-markov-final_bfstack.nii.gz")):
    print(" fully completed")
else:    
    stack_file = output3_dir+"foreground-markov-final_bfstack.nii.gz"
    bf = nighres.io.load_volume(blockface_img)
    stack_img = numpy.zeros(bf.header.get_data_shape())
        
    for idx,stain in enumerate(stains_nums):
        for num in stain:
            if int(num)>=first and int(num)<=last:
                histo = output3_dir+'foreground-markov-final-'+num+'_def-img.nii.gz'
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                
                if os.path.exists(histo):
                    stack_img[:,:,bfnum] = nighres.io.load_volume(histo).get_fdata()
    
    stack_nifti = nibabel.Nifti1Image(stack_img,affine=bf.affine,header=bf.header)
    nighres.io.save_volume(stack_file,stack_nifti)
     
     
# skullstrip, stack and intensity correct each set of stains
for idx,stain in enumerate(stains_nums):
    stack = []
    histo_stack = output3_dir+lb_stains[idx]+"stack.nii.gz"
    if (os.path.isfile(histo_stack)):
        print('Contrast '+str(idx)+' stacking: done')
    else:
        print('Contrast '+str(idx)+' stacking:')
        for num in stain:
            if int(num)>=first and int(num)<=last:
                histo = output3_dir+'foreground-markov-final-'+num+'_def-img.nii.gz'
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                mask = nighres.io.load_volume(blockface).get_fdata()>0
                
                stack.append(mask*nighres.io.load_volume(histo).get_fdata())
    
        img = numpy.stack(stack,axis=-1)
        header = nibabel.Nifti1Header()
        header.set_data_shape(img.shape)
    
        header.set_zooms((0.15,0.15,1.2))
        histo_nifti = nibabel.Nifti1Image(img,affine=None,header=header)
        nighres.io.save_volume(histo_stack,histo_nifti)

# stack corresponding masks
for idx,stain in enumerate(stains_nums):
    stack = []
    histo_stack = output3_dir+lb_stains[idx]+"fgstack.nii.gz"
    if (os.path.isfile(histo_stack)):
        print('Contrast '+str(idx)+' stacking: done')
    else:
        print('Contrast '+str(idx)+' stacking:')
        for num in stain:
            if int(num)>=first and int(num)<=last:
                histo = output3_dir+'foreground-markov6-'+num+'_ants-def1.nii.gz'
                
                if (int(num)<124): bfnum = int(num) + 16
                elif (int(num)<440): bfnum = int(num) + 17
                else: bfnum = int(num) + 16
                bfnum = str(bfnum).zfill(3)
                blockface = output1_dir+lb_blockface+bfnum+'-li-inv.nii.gz'
                mask = nighres.io.load_volume(blockface).get_fdata()>0
                
                stack.append(mask*nighres.io.load_volume(histo).get_fdata())
                #stack.append(nighres.io.load_volume(histo).get_fdata())
    
        img = numpy.stack(stack,axis=-1)
        header = nibabel.Nifti1Header()
        header.set_data_shape(img.shape)
    
        header.set_zooms((0.15,0.15,1.2))
        histo_nifti = nibabel.Nifti1Image(img,affine=None,header=header)
        nighres.io.save_volume(histo_stack,histo_nifti)

# intensity correct each set of stains 
for idx,stain in enumerate(stains_nums):
    stack = []
    histo_stack = output3_dir+lb_stains[idx]+"stack.nii.gz"
    histo_fgstack = output3_dir+lb_stains[idx]+"fgstack.nii.gz"
    nighres.microscopy.stack_intensity_regularisation(histo_stack, cutoff=75, mask=histo_fgstack, save_data=True)

# interpolate with nlm
for num,lb_stain in enumerate(lb_stains):
    histo_interp = output3_dir+lb_stain+'_foreground-markov_interp.nii.gz'
    if (first==0 and os.path.isfile(histo_interp)):
        print('Stain '+str(num)+' interpolation: done')
    else:
        print('Stain '+str(num)+' interpolation:')
        bf_img = nighres.io.load_volume(blockface_img).get_fdata()
        interp_img = numpy.zeros(bf_img.shape)
        
        histo_stack = nighres.io.load_volume(output3_dir+lb_stains[num]+"stack_sir-img.nii.gz").get_fdata()
    
        histo_images = stains_nums[num]
        for idx in range(bf_img.shape[2]):
            print("slice "+str(idx))
            # find closest stain slices
            prev = -1
            next = bf_img.shape[2]+1
            equal = -1
            idp = -1
            idn = -1
            ide = -1
            ndp = -1
            ndn = -1
            nde = -1
            for img,histo in enumerate(histo_images):
                
                # change the numbering to match
                if (int(histo)<124): bf_num = int(histo) + 16
                elif (int(histo)<440): bf_num = int(histo) + 17
                else: bf_num = int(histo) + 16

                if bf_num==idx: 
                    equal = bf_num
                    nde = img
                    ide = histo
                if bf_num>prev and bf_num<idx: 
                    prev = bf_num
                    ndp = img
                    idp = histo
                if bf_num<next and bf_num>idx: 
                    next = bf_num
                    ndn = img
                    idn = histo
            
            # check for borders
            if idn is histo_images[0]:
                next = bf_img.shape[2]+1
            if idp is histo_images[-1]:
                prev = -1
                    
            # if equal, use the stain directly
            if equal>-1:
                print("interpolate "+ide+" (in stack: "+str(nde))
                
                histo_curr = nibabel.Nifti1Image(histo_stack[:,:,nde],None,None)
                bf_curr = nibabel.Nifti1Image(bf_img[:,:,idx],None,None)
                
                nlm = nighres.microscopy.stack_intensity_mapping(image=bf_curr, 
                                        references=[bf_curr],
                                        mapped=[histo_curr],
                                        patch=5,search=5,
                                        save_data=True,overwrite=False,
                                        output_dir=output3_dir,
                                        file_name='interpolate_stain-'+lb_stains[num]+'_nlm-'+str(idx))
            
                interp_img[:,:,idx] = nighres.io.load_volume(nlm['result']).get_fdata()
            elif prev>-1 and next<bf_img.shape[2]+1:
                print("interpolate "+idp+" (in stack: "+str(ndp))
                print("and "+idn+" (in stack: "+str(ndn))
                
                histo_prev = nibabel.Nifti1Image(histo_stack[:,:,ndp],None,None)
                histo_next = nibabel.Nifti1Image(histo_stack[:,:,ndn],None,None)
                
                bf_curr = nibabel.Nifti1Image(bf_img[:,:,idx],None,None)
                bf_prev = nibabel.Nifti1Image(bf_img[:,:,prev],None,None)
                bf_next = nibabel.Nifti1Image(bf_img[:,:,next],None,None)
                
                wdist=[(next-idx)/(next-prev),(idx-prev)/(next-prev)]
                
                nlm = nighres.microscopy.stack_intensity_mapping(image=bf_curr, 
                                        references=[bf_prev,bf_next],
                                        mapped=[histo_prev,histo_next],
                                        weights=wdist,patch=5,search=5,
                                        save_data=True,overwrite=False,
                                        output_dir=output3_dir,
                                        file_name='interpolate_stain-'+lb_stains[num]+'_nlm-'+str(idx))
            
                interp_img[:,:,idx] = nighres.io.load_volume(nlm['result']).get_fdata()
            elif prev>-1:
                print("interpolate "+idp+" (in stack: "+str(ndp))
                                
                histo_prev = nibabel.Nifti1Image(histo_stack[:,:,ndp],None,None)
                
                bf_curr = nibabel.Nifti1Image(bf_img[:,:,idx],None,None)
                bf_prev = nibabel.Nifti1Image(bf_img[:,:,prev],None,None)
                
                nlm = nighres.microscopy.stack_intensity_mapping(image=bf_curr, 
                                        references=[bf_prev],
                                        mapped=[histo_prev],patch=5,search=5,
                                        save_data=True,overwrite=False,
                                        output_dir=output3_dir,
                                        file_name='interpolate_stain-'+lb_stains[num]+'_nlm-'+str(idx))
                
                interp_img[:,:,idx] = nlm['result'].get_fdata()
            elif next<bf_img.shape[2]+1:
                print("interpolate "+idn+" (in stack: "+str(ndn))
                
                histo_next = nibabel.Nifti1Image(histo_stack[:,:,ndn],None,None)
                
                bf_curr = nibabel.Nifti1Image(bf_img[:,:,idx],None,None)
                bf_next = nibabel.Nifti1Image(bf_img[:,:,next],None,None)
                
                nlm = nighres.microscopy.stack_intensity_mapping(image=bf_curr, 
                                        references=[bf_next],
                                        mapped=[histo_next],patch=5,search=5,
                                        save_data=True,overwrite=False,
                                        output_dir=output3_dir,
                                        file_name='interpolate_stain-'+lb_stains[num]+'_nlm-'+str(idx))
                
                interp_img[:,:,idx] = nlm['result'].get_fdata()
                
        histo_nifti = nibabel.Nifti1Image(interp_img,affine=nighres.io.load_volume(blockface_img).affine,
                                                     header=nighres.io.load_volume(blockface_img).header)
        nighres.io.save_volume(histo_interp,histo_nifti)

