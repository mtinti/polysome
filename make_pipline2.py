# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 11:52:11 2014

@author: mtinti-x
"""
import os

replace_list = eval(open('vars2.txt').read())
#template_file = open(os.path.join('templates','__test_template3.sh')).read()
#template_file = open(os.path.join('templates','__test_template_barcode.sh')).read()
#template_file = open(os.path.join('templates','trinity_template.sh')).read()
#template_file = open(os.path.join('templates','scallop.sh')).read()
#template_file = open(os.path.join('templates','transdecoder_template.sh')).read()
#template_file = open(os.path.join('templates','trinity_template.sh')).read()
#print(template_file)
#print(template_file)
template_file = open(os.path.join('templates','count_ht_template.sh')).read()
if __name__ == '__main__':
    run_all_content = ''
    for dictionary in replace_list:
        print(dictionary['g_version'],dictionary['base_fastq'],dictionary['experiment'])
        sh_script_name = 's_'+dictionary['experiment']+'.sh'
        run_all_content+='chmod +x '+sh_script_name+'\n'
        run_all_content+='qsub '+sh_script_name+'\n'
        run_all_content+='mv '+sh_script_name+' '+dictionary['experiment']+'\n'
        
        sh_script_content = template_file.replace('{experiment}',
        dictionary['experiment']).replace('{base_fastq}',
        dictionary['base_fastq']).replace('{g_version}',
        dictionary['g_version'])
        #.format(
            #g_version=dictionary['g_version'],
            #base_fastq=dictionary['base_fastq'],
            #experiment=dictionary['experiment'],
            #)

        open(sh_script_name,'w').write(sh_script_content)
        
        
  

    #run_all_content+='mv '+'run_all_'+dictionary['experiment']+'.sh'+' '+dictionary['experiment']+'\n'
    open('run_all.sh','w').write(run_all_content)
    
    

