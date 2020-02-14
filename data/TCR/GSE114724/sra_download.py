#!/usr/bin
#-*-coding:utf-8-*-

import os
import threading

def sra_download(sra_id,save_path):
    cmd_download = 'prefetch %s -o %s'%(sra_id,save_path)
    print(cmd_download)
    os.system(cmd_download)

def sra_download_fd(sra_id):
    command = "fastq-dump "+sra_id
    os.system(command)

def axel(url,save_file):
    command = "axel -n 10 -o %s %s"%(save_file,url)
    os.system(command)

def sra_to_fastq(sra_file,fastq_dir):
    command = "fastq-dump --gzip --split-3 -O %s %s"%(fastq_dir,sra_file)
    os.system(command)

def mkdir(dirname):
    command = 'mkdir -p '+dirname
    os.system(command)

class MyThread(threading.Thread):
    def __init__(self,sra_id,save_path):
        threading.Thread.__init__(self)
        self.sra_id = sra_id
        self.save_path = save_path
    def run(self):
        sra_download(self.sra_id,self.save_path)

def download_sra(sra_id_file,all_save_path):
    with open(sra_id_file,'r')as fin:
        lines = fin.readlines()
        tsk = []
        for line in lines[1::]:
            sra_id = line.strip().split(',')[0]
            c_dir = os.path.join(all_save_path,sra_id)
            mkdir(c_dir)
            save_sra_path = os.path.join(c_dir,'%s.sra'%(sra_id))
            #sra_download_fd(sra_id)
            if not os.path.exists(save_sra_path):
                sra_download(sra_id,save_sra_path)       
        #     try:
        #         tsk.append(MyThread(sra_id,save_sra_path))
        #     except:
        #         print('Error: unable to start thread')
        # for t in tsk:
        #     t.start()
        #     while True:
        #         if (len(threading.enumerate()) <= int(1)):
        #             break
        # for t in tsk:
        #     t.join()

def download_sra_to_fastq(sra_id_files,save_path):
    with open(sra_id_file,'r')as fin:
        lines = fin.readlines()
        tsk = []
        for line in lines[1:]:
            sra_id = line.strip().split(',')[0]
            c_dir = os.path.join(all_save_path,sra_id)
            save_sra_path = os.path.join(c_dir,'%s.sra'%(sra_id))
            sra_to_fastq(save_sra_path,c_dir)
            axel(url,save_file)

def download_fastq(ftp_file,save_path):
    with open(ftp_file,'r')as fin:
        lines = fin.readlines()
        tsk = []
        for line in lines[1:2]:
            line = line.strip().split('\t')
            sra_id = line[4]
            c_dir = os.path.join(save_path,sra_id)
            #save_sra_path = os.path.join(c_dir,'%s.sra'%(sra_id))
            mkdir(c_dir)
            fastq_url1 = line[9].split(';')[0]
            fastq_url2 = line[9].split(';')[1]
            axel(fastq_url1,c_dir)
            axel(fastq_url2,c_dir)


if __name__=='__main__':
    sra_id_files = '/100T/data/single_cell_sequence/GSE114724/SraRunTable.txt'
    save_path = "/100T/data/single_cell_sequence/TCR/GSE114724/download_fastq"
    #method1:forloop to download sra data,then sra to fastq
    # while 1==1:
    #     download_sra(sra_id_files,save_path)
    #download_fastq(sra_id_files,save_path)

    #method2:download fastq.gz in EUA according to the eua ftp
    ftp_file = "/100T/data/single_cell_sequence/TCR/GSE114724/PRJNA472381.txt"
    download_fastq(ftp_file,save_path)