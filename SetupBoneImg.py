# exec(open('C:/Users/ztoth/Desktop/SetupBoneImg.py').read())
import time
path = 'C:/Users/ztoth/Documents/MATLAB/PoreFigs/'
print("Use showPorosity('<filename.tif>') to display a .tif file")
print('The path to <filename.tif> is set to '+path)
print("Use iso() to rotate the current figure to an isometric view")
print("Use top() to rotate the current figure to an topside view")

def showPorosity(filenm):
    baselayer = importlayer(filename=path+filenm,importer='[ITK Importer]',mode='data',inputfiles_id='1')
    resizedlayer = transform(layerids=baselayer, origin=[0,0,0], spacing=[1,1,1], replace=True)
    time.sleep(0.1)
    maskedlayer = threshold(layerid=resizedlayer[0], lower_threshold='0',upper_threshold='0')
    rszlyr_vis = resizedlayer[0]+'::visible'
    msklyr_vis = maskedlayer+'::visible'
    set(stateid=rszlyr_vis,value='false')
    set(stateid=msklyr_vis,value='false')
    time.sleep(0.5)
    computeisosurface(layerid=maskedlayer,quality_factor='1',show='true')
    set(stateid='interface::python_console_visibility',value='false')
    
def iso():
    snap(viewerid='0')
    rotateview(stateid='viewer0::volume_view',axis='[0,1,0]',angle='45')
    rotateview(stateid='viewer0::volume_view',axis='[1,0,0]',angle='30')
    autoview(viewerid='0')
def top():
    snap(viewerid='0')
    rotateview(stateid='viewer0::volume_view',axis='[0,1,0]',angle='0')
    rotateview(stateid='viewer0::volume_view',axis='[1,0,0]',angle='0')
    autoview(viewerid='0')