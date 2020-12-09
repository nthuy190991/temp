import os
import glob

os.chdir('D:\\THNGU52\\')
args="--parse txyzic"
dir_data=os.getcwd()+'\\'

list_las=glob.glob('*.las')
for f in list_las:
    #f=list_las[0]

    print(f)
    filename_in=f
    filename_out=f[0:-4]+'.txt'
    print('Input:'+f+', Output:'+filename_out)
    
    os.chdir('C:\\Program Files\\QGIS')
    print("OSGeo4W.bat las2txt -i "+dir_data+filename_in+" -o "+dir_data+
              filename_out+" "+args)
    os.system("OSGeo4W.bat las2txt -i "+dir_data+filename_in+" -o "+dir_data+
              filename_out+" "+args)
    
