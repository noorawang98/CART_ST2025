#------wave equation---------
y(x, t) = A * sin(kx - ωt + φ)
rename.cart = c(
  '0'='CM_CD8',
  '3'='CM_CD8',
  '4'='',
  '6'='EFF')
#Final reference T cell
setwd("~/LJX/sc_in/")
immune.combind = readRDS("./rds/immune.combined_20230901.RDS")
immune.combined@active.ident = immune.combined$seurat_clusters
reference@active.ident =reference$seurat_clusters
celltype = c(
  '0'='CD8_Tem',
  '1'='CD8_Tem',
  '2'='CD8_Tem',
  '3'='CD8_Tpex',
  '4'='CD4_Th',
  '5'='CD8_Tem',
  '6'='CD4_Tn',
  '7'='CD8_Tem',
  '8'='CD8_Tem',
  '9'='CD4_Th',
  '10'='CD8_Tn/cm',
  '11'='CD8_Trm',
  '12'='CD4_Tn',
  "13"="NKT",
  '14'='CD8_Tex',
  '15'='CD8_Tem',
  '16'='CD8_Tem',
  '17'='Tcycling',
  '18'='γδT',
  '19'='γδT',
  '20'='CD8_Tex',
  "21"='Tumor',
  '22'='CD4_Tn/cm',
  "23"="Treg",
  '24'='Myeloid',
  '25'='CD4_Tn/cm',
  '26'='BC',
  '28'='γδT',
  '27'='Fibro',
  '29'='Myeloid',
  '30'='Tcycling',
  '31'='Fibro',
  '32'='γδT',
  '33'='EpC',
  '34'='Myeloid',
  '35'='Plasma',
  '36'='Endo')

celltype = factor(celltype,levels=c('Tumor',
                                    'CD8_Tn/cm',     
                                    'CD8_Trm',
                                    'CD8_Tem',
                                    'CD8_Tpex',
                                    'CD8_Tex',
                                    "CD4_Tn/cm",
                                    'CD4_Th',
                                    'Treg',
                                    'NKT',
                                    "γδT",
                                    'Tcycling',
                                    'Myeloid',
                                    'BC',
                                    'Plasma',
                                    'Fibro',
                                    'EpC',
                                    'Endo'
                                    
))

reference =RenameIdents(reference,celltype)
DefaultAssay(reference) <-'RNA'
reference = SCTransform(reference)
setwd("~/LJX/st/L13_cluster/")
st.integrated = readRDS("./rds/L13_st.integrated.RDS")
reference=readRDS('./rds/reference.RDS')

st.n9 = myCCAmap(obj.query = st.n9,obj.reference = reference)


cart = immune.combined@assays$SCT@counts["mCherry-CAR",] %>% as.numeric()

 egfp = immune.combined@assays$SCT@counts["EGFP",] %>% as.numeric()

 cd3 =immune.combined@assays$SCT@counts["Cd3d",] %>% as.numeric()
 
 myeloid = immune.combined@assays$SCT@counts["Itgam",]%>%as.numeric()
 # df =data.frame(x=immune.combined$spatial_x,y=immune.combined$spatial_y)
 df =data.frame(x=immune.combined$spatial_x,y=immune.combined$spatial_y,
                CD3=cd3,EGFP=egfp,CART=cart,Myeloid = myeloid,
                Group=immune.combined$Group)
 
 
 # "#E64B35FF" "#00A087FF"
 
 
 
 pdf("spatial_label_flor_group.pdf",width = 8,height = 4)
 ggplot(df,aes(x=x,y=y,fill=CART))+geom_point(shape=21,size=1)+scale_fill_viridis_c(option = "A")+
   facet_grid(.~Group)+cowplot::theme_map()
 ggplot(df,aes(x=x,y=y,fill=EGFP))+geom_point(shape=21,size=1)+scale_fill_viridis_c()+
   facet_grid(.~Group)+cowplot::theme_map()

 ggplot(df,aes(x=x,y=y,fill=CD3))+geom_point(shape=21,size=1)+
   scale_fill_gradientn(colours = rcolors::rcolors$NCV_manga)+
   facet_grid(.~Group)+cowplot::theme_map()
 ggplot(df,aes(x=x,y=y,fill=Myeloid))+geom_point(shape=21,size=1)+
   scale_fill_gradientn(colours = rcolors::rcolors$NCV_manga)+
   facet_grid(.~Group)+cowplot::theme_map()
 
 
 dev.off()

#------set installation and install dependent packages-------------
 options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=TRUE)
 
 # 镜像设置
 options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
 options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
 options("download.file.method"="libcurl")
 options("url.method"="libcurl")
 
 package=c("ggplot2", "BiocManager")
 for (pkg in package) {
    if (!requireNamespace(pkg, quietly = TRUE)){
       install.packages(pkg)
    }
 }
 
 config <- "https://bioconductor.org/config.yaml"
 readLines(config)
 # options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=T)
 bioc_package = c('Seurat','clusterProfiler','org.Mm.eg.db','TxDb.Hsapiens.UCSC.hg19.knownGene','org.Hs.eg.db')
 for (pkg in bioc_package) {
    if (!requireNamespace(pkg, quietly = TRUE)){
       BiocManager::install(pkg,ask = F,update = F)
    }
 }
 
#----CARD installation--------------------
 library(R.utils)
 chooseBioCmirror()
 # options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
 options(download.file.method = 'libcurl')
 options(url.method='libcurl')
 
 options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=TRUE)
 BiocManager::install ("TOAST") 
 
 devtools:: install_github ('xuranw/MuSiC')
 devtools::install_github('YingMa0107/CARD')

 
# # ------------hclcolors---------
#  $`Pastel 1`
#  [1] "#FFC5D0" "#FEC7C7" "#FAC9BE" "#F5CCB5" "#EFCFAF" "#E7D2AA" "#DED5A7" "#D4D8A7" "#C9DBAA" "#BFDDAF"
#  [11] "#B4DFB5" "#AAE1BD" "#A2E2C6" "#9CE2CF" "#99E2D8" "#9AE1E1" "#9EDFE9" "#A6DDF0" "#B0DAF6" "#BCD7FA"
#  [21] "#C8D3FC" "#D5D0FC" "#E0CCFA" "#EAC9F6" "#F2C7F1" "#F8C5EA" "#FDC4E2" "#FFC5D9"
#  
#  $`Dark 2`
#  [1] "#C87A8A" "#C57D7C" "#C1816E" "#BA8561" "#B28955" "#A88E4B" "#9D9246" "#909646" "#82994C" "#719C55"
#  [11] "#5F9F61" "#4AA16E" "#30A37C" "#06A389" "#00A396" "#00A2A3" "#00A0AE" "#319DB7" "#4E99BE" "#6794C4"
#  [21] "#7E8FC7" "#9189C7" "#A284C5" "#AF7FC0" "#BA7BB8" "#C179AF" "#C678A4" "#C87998"
#  
#  $`Dark 3`
#  [1] "#E16A86" "#DD706F" "#D57754" "#CC7E33" "#C18500" "#B38C00" "#A39200" "#909800" "#799D00" "#5CA200"
#  [11] "#2FA636" "#00A955" "#00AB6E" "#00AD85" "#00AD9A" "#00ACAD" "#00A9BE" "#00A5CC" "#009FD8" "#3097E1"
#  [21] "#6C8EE6" "#9183E6" "#AD79E3" "#C270DB" "#D169D0" "#DB64C1" "#E163B0" "#E3669C"
#  
#  $`Set 2`
#  [1] "#ED90A4" "#EA9493" "#E59882" "#DD9D71" "#D3A263" "#C8A857" "#BAAD50" "#ABB150" "#99B657" "#85B963"
#  [11] "#6FBC72" "#55BF82" "#33C192" "#00C1A2" "#00C1B2" "#00C0C1" "#00BDCE" "#32BAD9" "#5AB5E2" "#79AFE8"
#  [21] "#94A9EC" "#ACA2EC" "#BF9CE9" "#D096E4" "#DC91DB" "#E58ED0" "#EB8DC3" "#ED8EB4"
#  
#  $`Set 3`
#  [1] "#FFB3B5" "#FBB6A8" "#F4BA9B" "#ECBE90" "#E2C288" "#D6C783" "#C9CB82" "#BBCF85" "#ABD28C" "#9AD596"
#  [11] "#89D7A2" "#79D9AF" "#6BDABC" "#63D9C9" "#61D8D6" "#68D6E2" "#77D3EC" "#89D0F4" "#9DCBFA" "#B0C6FD"
#  [21] "#C3C1FE" "#D4BBFC" "#E2B7F8" "#EEB3F1" "#F7B0E8" "#FDAFDD" "#FFAFD0" "#FFB0C3"
#  
#  $Warm
#  [1] "#ABB065" "#AFAF64" "#B4AE64" "#B8AC65" "#BDAB66" "#C1A968" "#C4A86A" "#C8A66C" "#CBA56F" "#CFA373"
#  [11] "#D2A276" "#D5A07A" "#D79F7F" "#DA9D83" "#DC9C87" "#DE9B8C" "#DF9A91" "#E19895" "#E2979A" "#E3969F"
#  [21] "#E495A4" "#E495A8" "#E494AD" "#E494B2" "#E393B6" "#E393BA" "#E193BF" "#E093C3"
#  
#  $Cold
#  [1] "#ACA4E2" "#A5A6E2" "#9FA8E2" "#98AAE1" "#91ABE1" "#8AADE0" "#83AFDE" "#7BB1DD" "#74B2DB" "#6CB4D9"
#  [11] "#65B5D6" "#5DB7D3" "#56B8D1" "#4FB9CD" "#48BACA" "#43BBC6" "#3EBCC3" "#3ABCBF" "#38BDBB" "#37BDB6"
#  [21] "#39BEB2" "#3BBEAE" "#3FBEA9" "#44BEA5" "#49BEA0" "#4FBE9B" "#56BD97" "#5CBD92"
#  
#  $Harmonic
#  [1] "#C7A76C" "#C2A968" "#BCAB66" "#B5AD64" "#AEAF64" "#A7B166" "#9FB368" "#97B56C" "#8FB770" "#86B875"
#  [11] "#7DBA7B" "#74BB81" "#6ABC88" "#61BD8E" "#57BD95" "#4EBE9C" "#45BEA3" "#3EBEAA" "#39BEB1" "#37BDB7"
#  [21] "#39BCBE" "#3FBBC4" "#47BAC9" "#51B9CE" "#5BB7D3" "#67B5D7" "#72B3DA" "#7DB0DD"
#  
#  $Dynamic
#  [1] "#DB9D85" "#D4A179" "#CBA56E" "#C0AA67" "#B4AE64" "#A6B266" "#97B56C" "#86B875" "#75BB81" "#62BD8D"
#  [11] "#50BE9B" "#40BEA8" "#38BDB5" "#3CBCC1" "#4CB9CC" "#60B6D5" "#76B2DB" "#8BADE0" "#9FA8E2" "#B1A2E1"
#  [21] "#C19DDE" "#CD99D8" "#D796D0" "#DE94C6" "#E393BB" "#E494AE" "#E396A0" "#E09993"
#  
#  $Grays
#  [1] "#1B1B1B" "#242424" "#2D2D2D" "#363636" "#3F3F3F" "#484848" "#515151" "#5A5A5A" "#646464" "#6D6D6D"
#  [11] "#767676" "#808080" "#898989" "#929292" "#9B9B9B" "#A4A4A4" "#ADADAD" "#B5B5B5" "#BEBEBE" "#C6C6C6"
#  [21] "#CECECE" "#D6D6D6" "#DDDDDD" "#E4E4E4" "#EBEBEB" "#F1F1F1" "#F6F6F6" "#F9F9F9"
#  
#  $`Light Grays`
#  [1] "#474747" "#4E4E4E" "#565656" "#5E5E5E" "#656565" "#6D6D6D" "#747474" "#7B7B7B" "#838383" "#8A8A8A"
#  [11] "#919191" "#979797" "#9E9E9E" "#A5A5A5" "#ABABAB" "#B1B1B1" "#B7B7B7" "#BDBDBD" "#C2C2C2" "#C7C7C7"
#  [21] "#CCCCCC" "#D1D1D1" "#D5D5D5" "#D9D9D9" "#DCDCDC" "#DFDFDF" "#E1E1E1" "#E2E2E2"
#  
#  $`Blues 2`
#  [1] "#023FA5" "#2548A4" "#3650A5" "#4359A7" "#4F61A9" "#5969AC" "#6371AF" "#6C78B3" "#7580B6" "#7D87B9"
#  [11] "#868EBD" "#8E95C0" "#959CC3" "#9DA3C6" "#A4AAC9" "#ABB0CC" "#B2B6CF" "#B8BCD1" "#BEC1D4" "#C4C7D6"
#  [21] "#CACCD8" "#CFD0DA" "#D3D5DC" "#D8D8DE" "#DBDCE0" "#DFDFE1" "#E1E1E2" "#E2E2E2"
#  
#  $`Blues 3`
#  [1] "#00366C" "#003E75" "#00477F" "#005089" "#005995" "#0062A0" "#006BAC" "#0074B7" "#017DC3" "#2C86CA"
#  [11] "#438ECF" "#5596D5" "#659FDA" "#72A7DF" "#7FAFE4" "#8CB6E9" "#97BEEE" "#A2C5F2" "#ADCCF6" "#B7D3F9"
#  [21] "#C1DAFC" "#CBE0FF" "#D4E6FF" "#DCEBFF" "#E4F0FF" "#ECF4FF" "#F3F7FF" "#F9F9F9"
#  
#  $`Purples 2`
#  [1] "#3C2692" "#433291" "#4A3C92" "#514594" "#584D97" "#5F569B" "#675E9F" "#6E66A3" "#756EA7" "#7C76AB"
#  [11] "#847DAF" "#8B85B4" "#928DB8" "#9995BC" "#A09CC1" "#A7A4C5" "#AEABC9" "#B5B2CD" "#BCBAD1" "#C3C1D6"
#  [21] "#CAC8D9" "#D0CFDD" "#D6D5E1" "#DDDCE5" "#E2E2E8" "#E8E7EB" "#EDEDEE" "#F1F1F1"
#  
#  $`Purples 3`
#  [1] "#312271" "#3A2C7A" "#433684" "#4C3F8F" "#55489A" "#5E50A6" "#6759B2" "#7062BE" "#786BC6" "#8175CB"
#  [11] "#8A7FD0" "#9288D5" "#9B92DA" "#A39BDF" "#ABA4E4" "#B3ACE8" "#BBB5EC" "#C3BDF0" "#CAC5F3" "#D1CCF7"
#  [21] "#D8D4F9" "#DEDBFC" "#E4E1FD" "#EAE7FF" "#EFEDFF" "#F3F2FF" "#F7F6FD" "#F9F9F9"
#  
#  $`Reds 2`
#  [1] "#7F000D" "#86121E" "#8D222A" "#932F35" "#993A40" "#9F4549" "#A54F53" "#AA595D" "#B06366" "#B56C6F"
#  [11] "#BA7578" "#BF7E81" "#C3878A" "#C89092" "#CC999B" "#D0A1A3" "#D4AAAB" "#D7B2B3" "#DBBABB" "#DEC1C2"
#  [21] "#E1C9CA" "#E4D0D1" "#E7D7D7" "#E9DDDE" "#ECE3E4" "#EEE9E9" "#EFEDED" "#F1F1F1"
#  
#  $`Reds 3`
#  [1] "#69000C" "#770312" "#850718" "#940C1D" "#A21122" "#B11527" "#C0192C" "#CF1D31" "#DF2135" "#EE253A"
#  [11] "#F43A49" "#F94C58" "#FE5B65" "#FF6972" "#FF767E" "#FF8389" "#FF8F95" "#FF9BA0" "#FFA6AA" "#FFB1B5"
#  [21] "#FFBCBF" "#FFC6C8" "#FFD0D2" "#FFD9DA" "#FFE2E3" "#FFEAEA" "#FCF1F1" "#F6F6F6"
#  
#  $`Greens 2`
#  [1] "#006027" "#1A6632" "#2A6D3D" "#377346" "#42794F" "#4D8058" "#578661" "#608C6A" "#6A9272" "#73987B"
#  [11] "#7C9E83" "#85A48B" "#8DAA93" "#96B09B" "#9EB6A3" "#A6BCAA" "#AEC1B2" "#B6C7B9" "#BDCCC0" "#C4D1C7"
#  [21] "#CBD6CD" "#D2DBD4" "#D8DFDA" "#DFE4DF" "#E4E8E5" "#E9EBEA" "#EEEEEE" "#F1F1F1"
#  
#  $`Greens 3`
#  [1] "#004616" "#00501D" "#005B24" "#03652B" "#0A6F31" "#107937" "#15833D" "#198D43" "#1D9748" "#20A14E"
#  [11] "#3BA85B" "#4EB068" "#5EB775" "#6DBE81" "#7BC58C" "#88CC97" "#95D2A2" "#A1D8AC" "#ACDDB6" "#B7E2C0"
#  [21] "#C1E7C9" "#CBECD1" "#D4EFDA" "#DDF3E1" "#E6F6E9" "#EDF8EF" "#F4F9F5" "#F9F9F9"
#  
#  $Oslo
#  [1] "#FCFCFC" "#EFF2F8" "#E2E8F3" "#D5DDEF" "#C8D3EA" "#BBC9E5" "#AEBFE1" "#A1B6DC" "#94ACD7" "#86A2D3"
#  [11] "#7899CE" "#6990CA" "#5986C6" "#467DC2" "#3974B9" "#346BAB" "#30629D" "#2B598F" "#275182" "#224875"
#  [21] "#1E4068" "#19375B" "#152F4F" "#112842" "#0E2037" "#0A192B" "#07101D" "#040404"
#  
#  $`Purple-Blue`
#  [1] "#6B0077" "#6C1E7D" "#6D2F84" "#6F3C8B" "#714892" "#735398" "#755E9F" "#7768A5" "#7A71AC" "#7C7BB2"
#  [11] "#7F84B7" "#838DBD" "#8796C2" "#8B9FC7" "#90A7CC" "#95AFD0" "#9BB7D5" "#A1BED9" "#A7C6DD" "#AECCE0"
#  [21] "#B5D3E3" "#BDD9E7" "#C4DFE9" "#CCE4EC" "#D5E9EE" "#DDEDF0" "#E6F0F1" "#F1F1F1"
#  
#  $`Red-Purple`
#  [1] "#7D0112" "#861029" "#8F1C39" "#972648" "#9F3056" "#A73A63" "#AE436F" "#B54D7B" "#BB5687" "#C16092"
#  [11] "#C6699C" "#CB72A7" "#D07CB0" "#D485BA" "#D88EC3" "#DB97CB" "#DEA0D3" "#E1A9DA" "#E4B2E1" "#E6BAE7"
#  [21] "#E8C2ED" "#EACAF2" "#ECD2F6" "#EED9F9" "#EFDFFB" "#F1E6FC" "#F2EBFC" "#F2F0F6"
#  
#  $`Red-Blue`
#  [1] "#A93154" "#AC345E" "#B03767" "#B23A70" "#B53E79" "#B74381" "#B94789" "#BA4C90" "#BB5197" "#BC569E"
#  [11] "#BD5BA4" "#BD61AB" "#BD66B0" "#BC6CB6" "#BC71BB" "#BB77C0" "#BA7CC4" "#B982C9" "#B788CD" "#B68DD0"
#  [21] "#B592D4" "#B398D7" "#B29DDA" "#B1A2DC" "#B0A8DF" "#AFADE1" "#AEB2E3" "#AEB6E5"
#  
#  $`Purple-Orange`
#  [1] "#5B3794" "#663A96" "#703E98" "#7A429A" "#83469C" "#8B4B9E" "#9450A0" "#9B55A2" "#A35AA3" "#AA5FA5"
#  [11] "#B165A7" "#B86BA8" "#BF71AA" "#C577AC" "#CB7DAD" "#D184AF" "#D68AB1" "#DC91B3" "#E198B5" "#E69FB7"
#  [21] "#EAA6B9" "#EEADBC" "#F2B4BF" "#F6BCC2" "#F9C3C6" "#FBCBCA" "#FDD3D0" "#F8DCD9"
#  
#  $`Purple-Yellow`
#  [1] "#80146E" "#83247E" "#84328D" "#83409B" "#7E4EA7" "#765DAF" "#6C6AB5" "#6177BA" "#5483BE" "#478EC1"
#  [11] "#3A99C2" "#30A2C2" "#2BABC2" "#2FB4C1" "#3ABCBF" "#4AC3BD" "#5BCABB" "#6CD0BA" "#7ED5B8" "#8FDAB7"
#  [21] "#A0DFB7" "#B0E3B8" "#BFE7BA" "#CDEABE" "#DAEDC2" "#E5EFC8" "#EEF1CF" "#F5F2D8"
#  
#  $`Blue-Yellow`
#  [1] "#2D3184" "#24438A" "#1B5390" "#0F6297" "#01709E" "#007DA4" "#0089A9" "#0F95AD" "#21A0B1" "#32AAB5"
#  [11] "#42B3B8" "#51BCBA" "#60C4BC" "#6FCCBD" "#7ED3BF" "#8CD9C0" "#99DEC2" "#A6E3C3" "#B3E7C5" "#BEEAC7"
#  [21] "#C9EDC9" "#D2EFCC" "#DBF0CE" "#E2F1D2" "#E9F2D5" "#EEF2D9" "#F1F2DE" "#F3F1E4"
#  
#  $`Green-Yellow`
#  [1] "#006E37" "#18773F" "#307F48" "#418750" "#509058" "#5E9760" "#6B9F68" "#77A770" "#82AE78" "#8DB580"
#  [11] "#97BB87" "#A1C28F" "#AAC896" "#B3CD9D" "#BCD3A4" "#C4D8AB" "#CBDCB1" "#D2E1B7" "#D9E5BD" "#DFE9C3"
#  [21] "#E4ECC9" "#E9EFCE" "#EEF1D3" "#F1F3D8" "#F5F5DC" "#F7F6E1" "#F8F7E5" "#F9F7EA"
#  
#  $`Red-Yellow`
#  [1] "#7D0112" "#841815" "#8A2718" "#91331B" "#983E1F" "#9E4824" "#A55229" "#AB5C2E" "#B16533" "#B76F39"
#  [11] "#BD7840" "#C28146" "#C78A4D" "#CC9354" "#D19B5C" "#D5A463" "#D9AC6B" "#DDB573" "#E1BD7C" "#E4C484"
#  [21] "#E7CC8D" "#EAD396" "#ECDA9F" "#EEE0A9" "#F0E6B3" "#F2EBBD" "#F3F0CA" "#F2F1E4"
#  
#  $Heat
#  [1] "#8E063B" "#98223F" "#A23243" "#AB4147" "#B44D4A" "#BD594D" "#C56551" "#CD6F54" "#D47A56" "#DA8459"
#  [11] "#E08D5C" "#E5965E" "#E99F61" "#EDA763" "#F0AF66" "#F2B669" "#F4BD6C" "#F5C36E" "#F6C971" "#F6CF75"
#  [21] "#F5D478" "#F4D97B" "#F2DD7F" "#F0E083" "#EDE388" "#EAE58E" "#E6E795" "#E2E6BD"
#  
#  $`Heat 2`
#  [1] "#D33F6A" "#D54766" "#D84F62" "#DA565E" "#DC5D59" "#DE6455" "#E06B50" "#E2714B" "#E37846" "#E57E41"
#  [11] "#E6853C" "#E78B37" "#E89132" "#E9972E" "#E99D2A" "#EAA428" "#EAAA28" "#EAB02A" "#E9B62D" "#E9BC33"
#  [21] "#E8C23A" "#E8C842" "#E7CE4B" "#E6D356" "#E5D961" "#E3DF6E" "#E2E47E" "#E2E6BD"
#  
#  $Terrain
#  [1] "#26A63A" "#43A831" "#57AA29" "#68AC20" "#77AE16" "#85B00D" "#92B106" "#9EB307" "#A9B40F" "#B4B61A"
#  [11] "#BFB726" "#C9B831" "#D3BA3D" "#DDBB49" "#E6BC54" "#EFBD60" "#F7BF6C" "#FFC078" "#FFC183" "#FFC38F"
#  [21] "#FFC49B" "#FFC6A7" "#FFC8B2" "#FFCBBE" "#FFCDC9" "#FFD1D5" "#FFD5E0" "#F1F1F1"
#  
#  $`Terrain 2`
#  [1] "#027C1E" "#2A8121" "#3F8626" "#4F8B2B" "#5D9030" "#6A9537" "#769A3D" "#829E44" "#8DA24C" "#97A753"
#  [11] "#A1AB5B" "#AAAF63" "#B3B36B" "#BBB673" "#C3BA7C" "#CBBE84" "#D2C18C" "#D8C494" "#DEC79D" "#E4CBA5"
#  [21] "#E8CEAD" "#ECD1B4" "#EFD3BC" "#F1D6C3" "#F2D9CA" "#F2DBD1" "#EFDED8" "#E2E2E2"
#  
#  $Viridis
#  [1] "#4B0055" "#4A0C5E" "#471E67" "#422C70" "#3A3878" "#2E4480" "#185086" "#005B8C" "#006591" "#007094"
#  [11] "#007A97" "#008498" "#008E98" "#009796" "#00A094" "#00A890" "#00B08B" "#00B784" "#00BE7D" "#05C574"
#  [21] "#4ACB6A" "#6CD05E" "#89D552" "#A3DA44" "#BBDD38" "#D2E02E" "#E9E22C" "#FDE333"
#  
#  $Plasma
#  [1] "#001889" "#1F0F89" "#42048A" "#58008B" "#6A008C" "#79008D" "#87008D" "#94008C" "#A0038B" "#AB1488"
#  [11] "#B52285" "#BE2F81" "#C73B7B" "#CE4875" "#D5546D" "#DB6063" "#E06C58" "#E5794B" "#E8853A" "#EB9223"
#  [21] "#EC9F00" "#EDAC00" "#ECB900" "#EBC600" "#E8D400" "#E4E200" "#DFF019" "#DAFF47"
#  
#  $Inferno
#  [1] "#040404" "#0D0C22" "#1A1032" "#27123E" "#36134A" "#451354" "#55125D" "#651165" "#75106B" "#851170"
#  [11] "#941473" "#A31A75" "#B12275" "#BE2C72" "#CB376E" "#D74366" "#E2505B" "#EC5E49" "#F36E35" "#F47F29"
#  [21] "#F69022" "#F7A026" "#F7B033" "#F8BF46" "#F9CE5A" "#FADE70" "#FCED86" "#FFFE9E"
#  
#  $Rocket
#  [1] "#070707" "#1B0921" "#2A072F" "#370439" "#440142" "#51004A" "#5F0052" "#6C0058" "#7A005D" "#870061"
#  [11] "#950064" "#A20065" "#AF0065" "#BC0464" "#C91660" "#D6245A" "#E23250" "#E6484A" "#E85D48" "#EA6F4B"
#  [21] "#EC7F51" "#EE8F5C" "#EF9F6A" "#F0AE7C" "#F2BE91" "#F4CEA8" "#F7DFC4" "#FDF5EB"
#  
#  $Mako
#  [1] "#070707" "#1B0D19" "#261225" "#2E1831" "#361E3D" "#3C2649" "#402E56" "#433763" "#43406F" "#404B7B"
#  [11] "#395687" "#2B6192" "#036D9B" "#0079A4" "#0086AC" "#0092B2" "#009DB4" "#00A6B3" "#00AFB2" "#00B8B3"
#  [21] "#3FC0B4" "#60C8B6" "#7BD0B9" "#92D8BE" "#A8DFC5" "#BCE7CD" "#CEEFD6" "#E0F7E1"
#  
#  $`Dark Mint`
#  [1] "#0E3F5C" "#114660" "#144C65" "#185369" "#1C5A6E" "#216172" "#266877" "#2C6F7B" "#31767F" "#387D84"
#  [11] "#3E8488" "#458B8D" "#4C9291" "#539995" "#5BA09A" "#62A79E" "#6AAEA2" "#73B5A6" "#7BBCAB" "#84C3AF"
#  [21] "#8DCAB3" "#96D1B8" "#A0D8BC" "#A9DFC1" "#B3E6C5" "#BDEDCA" "#C7F4CF" "#D1FBD4"
#  
#  $Mint
#  [1] "#005D67" "#00626B" "#00686E" "#006D71" "#007275" "#007878" "#0A7D7C" "#1F8380" "#2D8884" "#398E88"
#  [11] "#43948B" "#4D998F" "#569F94" "#5FA598" "#68AA9C" "#71B0A1" "#7AB6A5" "#82BBAA" "#8BC1AF" "#94C7B4"
#  [21] "#9CCCB9" "#A5D2BE" "#AED8C3" "#B7DDC9" "#C0E3CF" "#C9E9D5" "#D3EEDC" "#E0F2E6"
#  
#  $BluGrn
#  [1] "#14505C" "#175661" "#195C65" "#1D6269" "#20686E" "#236E72" "#277576" "#2B7B7A" "#2F817D" "#348781"
#  [11] "#398D84" "#3E9488" "#449A8B" "#49A08E" "#50A690" "#56AC93" "#5DB295" "#63B897" "#6BBE99" "#72C49B"
#  [21] "#7EC89E" "#8ACDA2" "#95D2A6" "#A0D6AA" "#AADBAE" "#B4DFB3" "#BDE2B8" "#C7E5BE"
#  
#  $Teal
#  [1] "#2A5676" "#2D5C7C" "#2F6381" "#326986" "#346F8C" "#377591" "#397C96" "#3C829B" "#4288A0" "#498EA4"
#  [11] "#4F94A8" "#569AAC" "#5DA0B0" "#64A6B4" "#6BACB8" "#72B2BC" "#79B7C0" "#80BDC4" "#88C3C8" "#8FC8CC"
#  [21] "#97CED0" "#9ED3D4" "#A6D8D8" "#AEDDDC" "#B6E2E0" "#BEE7E3" "#C6EBE7" "#D2EEEA"
#  
#  $TealGrn
#  [1] "#177F97" "#11859A" "#098A9D" "#0090A0" "#0095A2" "#009AA4" "#00A0A7" "#00A5A8" "#00AAAA" "#00AFAB"
#  [11] "#09B5AD" "#15BAAE" "#20BEAE" "#29C3AF" "#33C8AF" "#3CCDAF" "#45D1AF" "#4ED5AF" "#57DAAE" "#61DDAD"
#  [21] "#6EE1AD" "#7BE4AE" "#86E7AE" "#91EAAE" "#9CECAF" "#A6EEB0" "#AFF0B1" "#B7F1B2"
#  
#  $Emrld
#  [1] "#0B4151" "#074756" "#024E5B" "#005560" "#005C64" "#006269" "#00696D" "#007071" "#007775" "#007E78"
#  [11] "#09857B" "#168C7E" "#229380" "#2D9A83" "#38A085" "#43A787" "#4EAE89" "#59B58B" "#65BB8C" "#71C28E"
#  [21] "#7DC890" "#89CE91" "#95D593" "#A1DB96" "#AEE198" "#BBE79B" "#C7ED9F" "#D4F3A3"
#  
#  $BluYl
#  [1] "#2E4F79" "#29577E" "#235E81" "#1C6685" "#156D89" "#0C758C" "#057C8E" "#028391" "#058A93" "#109195"
#  [11] "#1B9897" "#279F98" "#33A69A" "#3FAD9B" "#4CB39C" "#58BA9D" "#65C09D" "#72C69E" "#7FCD9F" "#8CD3A0"
#  [21] "#99D9A1" "#A6DEA2" "#B4E4A4" "#C2EAA6" "#CFEFA8" "#DDF5AB" "#EBFAAE" "#F9FFAF"
#  
#  $ag_GrnYl
#  [1] "#255668" "#1F5C6C" "#186270" "#0C6874" "#006E77" "#007479" "#007A7C" "#00807D" "#00867F" "#008C80"
#  [11] "#009280" "#009880" "#009E80" "#07A47F" "#22A97E" "#34AF7C" "#44B57A" "#53BA77" "#61C074" "#70C571"
#  [21] "#7FCA6D" "#8ECF69" "#9DD566" "#ACDA62" "#BBDF5F" "#CBE45D" "#DBE95B" "#EDEF5C"
#  
#  $Peach
#  [1] "#EA4C3B" "#EB543D" "#EC5C40" "#ED6343" "#EE6946" "#EF704A" "#F0764E" "#F17C52" "#F28257" "#F3885B"
#  [11] "#F48E60" "#F59366" "#F5996B" "#F69E70" "#F6A376" "#F7A87C" "#F7AD82" "#F8B288" "#F8B78E" "#F8BC94"
#  [21] "#F9C19A" "#F9C5A0" "#F9CAA6" "#FACEAC" "#FAD2B3" "#FAD6B8" "#FADABE" "#FADDC3"
#  
#  $PinkYl
#  [1] "#E24C80" "#E5547D" "#E75C7A" "#E86478" "#EA6B76" "#EC7274" "#ED7972" "#EF8071" "#F0866F" "#F28D6F"
#  [11] "#F3936F" "#F4996F" "#F5A070" "#F5A671" "#F6AC73" "#F7B275" "#F7B878" "#F8BE7C" "#F9C480" "#F9C984"
#  [21] "#F9CF89" "#FAD58F" "#FADA94" "#FBE09A" "#FBE6A1" "#FCEBA7" "#FDF0AE" "#FDF6B5"
#  
#  $Burg
#  [1] "#67223F" "#6F2644" "#762A48" "#7E2E4D" "#863251" "#8E3655" "#963B5A" "#9E3F5E" "#A64362" "#AE4867"
#  [11] "#B64C6B" "#BE516F" "#C55673" "#C95E78" "#CD667D" "#D16E82" "#D57587" "#D97D8D" "#DD8492" "#E18C98"
#  [21] "#E5939D" "#E99AA3" "#EDA1A9" "#F1A9AF" "#F6B0B5" "#FAB7BB" "#FEBEC1" "#FFC5C7"
#  
#  $BurgYl
#  [1] "#772C4B" "#7E314E" "#853650" "#8C3B53" "#934155" "#9A4657" "#A14B59" "#A8515A" "#AF565C" "#B75C5D"
#  [11] "#BE615E" "#C5675E" "#CC6D5E" "#D3735D" "#D9795D" "#DE815F" "#E18963" "#E49168" "#E7996E" "#E9A174"
#  [21] "#ECA87B" "#EEB083" "#F1B88B" "#F3C093" "#F5C79C" "#F7CFA7" "#F8D7B2" "#F8DFC1"
#  
#  $RedOr
#  [1] "#B13F63" "#B64463" "#BB4863" "#C04C63" "#C55063" "#CA5562" "#CF5961" "#D45E60" "#D8635F" "#DA695F"
#  [11] "#DC7060" "#DF7661" "#E07C63" "#E28265" "#E48867" "#E68E6A" "#E8936D" "#E99970" "#EB9F74" "#ECA578"
#  [21] "#EEAA7C" "#EFB081" "#F1B586" "#F2BB8C" "#F3C192" "#F4C698" "#F5CCA0" "#F5D1A8"
#  
#  $OrYel
#  [1] "#EF4868" "#F34C62" "#F5525D" "#F55A59" "#F56155" "#F56852" "#F56F4F" "#F5764C" "#F57C4A" "#F58249"
#  [11] "#F58848" "#F48D48" "#F49349" "#F3984B" "#F39D4D" "#F2A251" "#F2A855" "#F1AC59" "#F0B15E" "#F0B664"
#  [21] "#EFBB6A" "#EEBF70" "#EEC476" "#EDC87D" "#EDCD84" "#ECD18B" "#ECD592" "#ECD999"
#  
#  $Purp
#  [1] "#645A9F" "#6A5EA4" "#7162AA" "#7766AF" "#7E69B5" "#846DBA" "#8A72BD" "#9077C0" "#967BC4" "#9C80C7"
#  [11] "#A285CB" "#A88ACE" "#AE8FD1" "#B393D5" "#B998D8" "#BF9EDB" "#C4A3DF" "#C9A8E2" "#CFADE5" "#D4B2E8"
#  [21] "#D9B8EC" "#DEBDEF" "#E3C3F2" "#E8C8F5" "#EDCEF7" "#F1D4FA" "#F5DAFC" "#F6E2FB"
#  
#  $PurpOr
#  [1] "#583C88" "#633E8D" "#6E4192" "#784496" "#82479B" "#8C4A9F" "#964DA3" "#A051A6" "#A855A9" "#AF5BAA"
#  [11] "#B661AB" "#BD67AC" "#C36EAD" "#C974AE" "#CF7BAF" "#D482B0" "#DA89B2" "#DE90B3" "#E397B5" "#E79EB7"
#  [21] "#ECA5B9" "#EFACBC" "#F3B4BF" "#F6BCC2" "#F8C3C6" "#FACBCB" "#FCD3D0" "#F8DCD9"
#  
#  $Sunset
#  [1] "#704D9E" "#7D4EA1" "#8A4FA4" "#9650A6" "#A152A8" "#AB55A9" "#B557A9" "#BE5BA9" "#C75EA8" "#CF63A6"
#  [11] "#D768A4" "#DE6DA1" "#E4739D" "#EA7999" "#EF8095" "#F28891" "#F4908D" "#F69889" "#F7A086" "#F8A884"
#  [21] "#F9B082" "#F9B881" "#F9C082" "#F8C884" "#F7D087" "#F6D88C" "#F5DF92" "#F3E79A"
#  
#  $Magenta
#  [1] "#6D1C68" "#74226C" "#7B2971" "#822E76" "#89347A" "#913A7F" "#983F84" "#9F4488" "#A7498D" "#AF4F91"
#  [11] "#B65496" "#BE599A" "#C3619D" "#C868A1" "#CD6FA4" "#D276A7" "#D77DAA" "#DB84AE" "#E08BB1" "#E492B5"
#  [21] "#E899B8" "#ECA0BC" "#EFA7BF" "#F2AEC3" "#F5B5C7" "#F7BBCB" "#F7C2CE" "#F3CAD2"
#  
#  $SunsetDark
#  [1] "#7D1D67" "#881E6B" "#921E6F" "#9C2073" "#A72275" "#B12477" "#BB2879" "#C42C7A" "#CE3079" "#D73679"
#  [11] "#E03B77" "#E94274" "#F14871" "#F5546E" "#F7606D" "#F96C6C" "#FB776C" "#FD816D" "#FE8B6F" "#FF9471"
#  [21] "#FF9E75" "#FFA779" "#FFAF7F" "#FFB884" "#FFC18B" "#FFC992" "#FFD198" "#FFD99F"
#  
#  $ag_Sunset
#  [1] "#4B1D91" "#5D1994" "#6D1597" "#7B1299" "#88109B" "#940F9C" "#A0119D" "#AB159C" "#B51A9B" "#BF219A"
#  [11] "#C92997" "#D13194" "#DA3A90" "#E1428A" "#E94C84" "#ED577C" "#EF6475" "#F1706E" "#F27B68" "#F38663"
#  [21] "#F3905F" "#F39B5D" "#F2A55E" "#F1AF62" "#EFB869" "#EDC274" "#EBCB82" "#E7D39A"
#  
#  $BrwnYl
#  [1] "#541C3A" "#5C223F" "#642843" "#6D2E47" "#76344A" "#7E3A4E" "#874151" "#904754" "#984E56" "#A15558"
#  [11] "#A95C5A" "#B1635B" "#B76B5D" "#BB7462" "#BF7D67" "#C3856C" "#C78E72" "#CB9779" "#CE9F80" "#D2A787"
#  [21] "#D5AF8F" "#D8B797" "#DCBFA0" "#DFC7A8" "#E2CEB1" "#E5D5BA" "#E9DCC2" "#EBE2CA"
#  
#  $YlOrRd
#  [1] "#7D0025" "#8C0020" "#9B0413" "#A90D00" "#B71800" "#C42300" "#D12D00" "#DD3800" "#E84300" "#EB5500"
#  [11] "#ED6500" "#EF7300" "#F18000" "#F38D00" "#F49800" "#F5A400" "#F6AE00" "#F6B836" "#F7C252" "#F8CB68"
#  [21] "#F9D47A" "#FADC8A" "#FBE399" "#FCEAA6" "#FEF1B2" "#FFF7BC" "#FFFCC3" "#FFFFC8"
#  
#  $YlOrBr
#  [1] "#682714" "#772D0B" "#863500" "#953C00" "#A34400" "#B04D00" "#BD5500" "#CA5E00" "#D66700" "#E17100"
#  [11] "#E97C00" "#EC8800" "#EF9300" "#F19F00" "#F3A900" "#F5B300" "#F6BD2C" "#F8C650" "#F8CF69" "#F9D77E"
#  [21] "#FADE90" "#FBE5A1" "#FCEBB0" "#FCF1BE" "#FDF5CA" "#FEF9D5" "#FEFCDE" "#FEFEE3"
#  
#  $OrRd
#  [1] "#88002D" "#940030" "#A10032" "#AE0732" "#BB1A32" "#C7272F" "#D4332A" "#E03E21" "#E74E26" "#EC5D2F"
#  [11] "#F16A38" "#F57742" "#F9834C" "#FC8E57" "#FF9961" "#FFA46C" "#FFAE77" "#FFB782" "#FFC08E" "#FFC999"
#  [21] "#FFD1A4" "#FFD8AF" "#FFDFBA" "#FFE5C5" "#FFEBD0" "#FFEFDA" "#FFF3E3" "#FDF5EB"
#  
#  $Oranges
#  [1] "#802A07" "#8C2F00" "#993400" "#A63A00" "#B23F00" "#BE4500" "#CA4B00" "#D65200" "#E25800" "#E96200"
#  [11] "#ED6D00" "#F07800" "#F48300" "#F68D00" "#F99730" "#FBA04A" "#FCAA5E" "#FEB370" "#FFBB80" "#FFC48F"
#  [21] "#FFCB9E" "#FFD3AB" "#FFDAB8" "#FFE1C5" "#FFE7D0" "#FFEDDB" "#FFF2E5" "#FEF5EC"
#  
#  $YlGn
#  [1] "#004533" "#00503A" "#005B40" "#006646" "#02714C" "#117C50" "#1E8654" "#2A9058" "#359A5B" "#40A45D"
#  [11] "#4CAE5E" "#57B75F" "#6BBE67" "#7DC571" "#8DCB7B" "#9CD185" "#AAD78F" "#B6DD99" "#C2E2A2" "#CCE7AC"
#  [21] "#D6EBB5" "#DFF0BE" "#E7F3C6" "#EEF7CD" "#F3F9D5" "#F8FCDB" "#FCFDE0" "#FEFEE3"
#  
#  $YlGnBu
#  [1] "#26185F" "#1D286E" "#00377D" "#00468B" "#005598" "#0064A4" "#0073AD" "#0081B2" "#008BB0" "#0095AF"
#  [11] "#009FAE" "#00A8AE" "#00B1AE" "#00B9AF" "#37C1B0" "#59C8B2" "#72CFB5" "#87D6B8" "#9ADCBB" "#AAE2C0"
#  [21] "#B9E8C4" "#C6EDC8" "#D2F2CD" "#DDF6D1" "#E6F9D5" "#EFFCD9" "#F6FEDC" "#FCFFDD"
#  
#  $Reds
#  [1] "#6D0026" "#7B002A" "#89002E" "#970031" "#A50034" "#B30A36" "#C11437" "#CF1C36" "#DD2434" "#EB2C31"
#  [11] "#F23E38" "#F75043" "#FB604F" "#FF6E5A" "#FF7B65" "#FF8870" "#FF947C" "#FF9F87" "#FFAA93" "#FFB59E"
#  [21] "#FFBFAA" "#FFC8B5" "#FFD2C1" "#FFDACC" "#FFE2D7" "#FFEAE1" "#FFF0EB" "#FCF5F2"
#  
#  $RdPu
#  [1] "#490062" "#550069" "#620071" "#710479" "#800B82" "#8F128A" "#9E1892" "#AE1E99" "#BD249F" "#CC2AA5"
#  [11] "#DA33A9" "#E048A8" "#E659A8" "#EB69A9" "#EF77AA" "#F484AB" "#F790AE" "#FB9CB1" "#FDA8B4" "#FFB3B9"
#  [21] "#FFBDBE" "#FFC7C4" "#FFD1CB" "#FFDAD2" "#FFE2DA" "#FFEAE2" "#FFF0EB" "#FCF5F2"
#  
#  $PuRd
#  [1] "#611300" "#710E00" "#82040A" "#930021" "#A40034" "#B60046" "#C7005A" "#D7006E" "#DB007F" "#DD008F"
#  [11] "#DE209C" "#DE3CA8" "#DE51B3" "#DC62BC" "#DB72C4" "#D981CB" "#D78ED1" "#D69BD6" "#D5A6DB" "#D6B2DF"
#  [21] "#D7BCE2" "#D9C6E6" "#DCCFE9" "#E0D8EC" "#E5E1F0" "#EBE8F3" "#F0F0F8" "#F6F6FC"
#  
#  $Purples
#  [1] "#3D1778" "#452380" "#4D2D8A" "#553695" "#5D3FA0" "#644BA3" "#6B56A7" "#7360AB" "#7B6AB0" "#8273B5"
#  [11] "#8A7DBA" "#9286BF" "#9A8FC4" "#A298CA" "#A9A1CF" "#B1AAD4" "#B9B3D8" "#C0BBDD" "#C8C3E2" "#CFCBE6"
#  [21] "#D6D3EB" "#DDDAEF" "#E3E1F3" "#E9E8F7" "#EFEEFA" "#F5F3FD" "#F9F8FF" "#FCFBFF"
#  
#  $PuBuGn
#  [1] "#004533" "#004F3D" "#005848" "#006155" "#006A62" "#007370" "#007C7F" "#00848F" "#008CA0" "#0091AD"
#  [11] "#0096B6" "#139ABF" "#4C9FC6" "#6AA3CD" "#81A8D2" "#94AED7" "#A5B3DC" "#B3B9DF" "#C0BFE3" "#CBC5E6"
#  [21] "#D5CCE8" "#DED3EB" "#E5DAEE" "#ECE0F1" "#F2E7F4" "#F7EDF7" "#FBF3FA" "#FFF7FD"
#  
#  $PuBu
#  [1] "#0E3F5C" "#0A4769" "#074E75" "#055682" "#055E8F" "#07659B" "#0C6DA8" "#1275B5" "#2F7CBC" "#4782BF"
#  [11] "#5A89C2" "#6990C5" "#7797C9" "#839ECC" "#8FA6D0" "#9AADD3" "#A4B4D7" "#AEBBDB" "#B7C2DE" "#C0C9E2"
#  [21] "#C9D0E6" "#D1D7E9" "#D9DDED" "#E1E4F1" "#E8EAF5" "#EEF0F8" "#F4F5FC" "#F8F9FF"
#  
#  $Greens
#  [1] "#004616" "#00501D" "#075A23" "#126429" "#1B6E2E" "#237833" "#2B8238" "#328C3C" "#399540" "#419F44"
#  [11] "#4EA74D" "#5EAF5B" "#6DB667" "#7ABD74" "#87C47F" "#93CA8B" "#9FD196" "#AAD7A1" "#B4DCAB" "#BEE2B5"
#  [21] "#C7E7BE" "#D0EBC8" "#D8EFD1" "#E0F3D9" "#E7F6E1" "#EDF8E8" "#F2FAEF" "#F6FBF4"
#  
#  $BuGn
#  [1] "#1B4414" "#1C4F1C" "#1C5A25" "#19662E" "#127137" "#007C41" "#00874B" "#009255" "#009C62" "#04A472"
#  [11] "#2EAC80" "#45B48E" "#58BB9A" "#69C2A6" "#78C9B1" "#87CFBC" "#94D5C5" "#A1DACE" "#ADE0D6" "#B9E4DD"
#  [21] "#C4E9E4" "#CDEDEA" "#D6F0EF" "#DEF4F3" "#E6F6F6" "#EBF9F9" "#F0FAFB" "#F2FBFC"
#  
#  $GnBu
#  [1] "#2F327D" "#23418A" "#005096" "#005EA1" "#006DAB" "#007AB4" "#0087BA" "#0092BA" "#009DBB" "#00A6BB"
#  [11] "#00AFBC" "#00B8BC" "#1AC0BD" "#45C7BD" "#60CEBF" "#76D4C0" "#89DAC2" "#9ADFC5" "#A9E3C8" "#B7E8CB"
#  [21] "#C4EBCF" "#CFEED3" "#D9F1D8" "#E1F4DD" "#E8F6E1" "#EEF7E5" "#F2F8E8" "#F5F8EA"
#  
#  $BuPu
#  [1] "#540046" "#5B0052" "#630D5F" "#6B1B6C" "#72277A" "#793388" "#7F3E96" "#834AA3" "#8459AA" "#8566B1"
#  [11] "#8672B7" "#887EBD" "#8B89C2" "#8E94C7" "#939FCC" "#98A9D1" "#9EB2D6" "#A4BBDA" "#ABC4DE" "#B3CCE2"
#  [21] "#BBD4E6" "#C4DBE9" "#CDE2ED" "#D5E9F0" "#DEEFF4" "#E6F4F7" "#EDF8FA" "#F2FBFC"
#  
#  $Blues
#  [1] "#273871" "#2B407B" "#2E4986" "#305292" "#315B9D" "#3264A8" "#316DB3" "#3976BA" "#457EBE" "#5087C1"
#  [11] "#5B8FC5" "#6697C9" "#709FCD" "#7AA7D1" "#85AFD5" "#8FB7D9" "#99BEDD" "#A2C5E0" "#ACCCE4" "#B6D3E8"
#  [21] "#BFD9EB" "#C8DFEE" "#D1E5F2" "#D9EAF5" "#E1F0F8" "#E9F4FA" "#EFF8FC" "#F4FAFE"
#  
#  $Lajolla
#  [1] "#FCFFC9" "#F8F8BA" "#F4EFAC" "#F1E59D" "#EEDC8F" "#ECD280" "#EAC872" "#E7BE63" "#E5B455" "#E3A946"
#  [11] "#E09E36" "#DE9326" "#DB8712" "#D87B00" "#D27004" "#CA6623" "#C05C30" "#B65338" "#AB4A3D" "#A04140"
#  [21] "#943841" "#873040" "#7A283E" "#6C203A" "#5D1935" "#4E122E" "#3C0D25" "#1D0B14"
#  
#  $Turku
#  [1] "#FFE9EA" "#FFDFDE" "#FFD6D1" "#FFCDC3" "#FCC4B5" "#F7BBA6" "#F0B395" "#E7AC86" "#D9A67F" "#CBA079"
#  [11] "#BE9A73" "#B1936D" "#A48D68" "#988664" "#8C7E5F" "#80775B" "#756F56" "#6B6751" "#605F4C" "#575747"
#  [21] "#4D4F42" "#44463C" "#3B3D36" "#32342F" "#292B27" "#20211F" "#171716" "#040404"
#  
#  $Hawaii
#  [1] "#8B0069" "#910C63" "#971E5C" "#9C2C55" "#A0384C" "#A34341" "#A64D35" "#A85724" "#A96204" "#A96C00"
#  [11] "#A77600" "#A58000" "#A18A00" "#9C9300" "#979D00" "#90A61C" "#88AF39" "#7FB850" "#75C165" "#6AC979"
#  [21] "#60D18D" "#57D8A0" "#51DFB2" "#51E5C4" "#59EBD5" "#68F0E5" "#7EF4F3" "#B0F4FA"
#  
#  $Batlow
#  [1] "#201158" "#01235A" "#002F5C" "#003A5E" "#004560" "#004E60" "#005760" "#00605D" "#00685A" "#007054"
#  [11] "#00774D" "#007D44" "#228339" "#48892A" "#648D16" "#7D9200" "#949500" "#AB9800" "#C19A1B" "#D49C41"
#  [21] "#E49E64" "#F2A181" "#FDA49B" "#FFAAB4" "#FFB0C9" "#FFB8DC" "#FFC2EC" "#FFCEF4"
#  
#  $`Blue-Red`
#  [1] "#023FA5" "#3650A5" "#4F61A9" "#6371AF" "#7580B6" "#868EBD" "#959CC3" "#A4AAC9" "#B2B6CF" "#BEC1D4"
#  [11] "#CACCD8" "#D3D5DC" "#DBDCE0" "#E1E1E2" "#E2E1E1" "#E0DBDC" "#DED2D4" "#DAC8CB" "#D6BCC0" "#D1AEB4"
#  [21] "#CCA0A8" "#C6909A" "#BE7F8C" "#B66E7D" "#AE5A6D" "#A4465D" "#9A2E4C" "#8E063B"
#  
#  $`Blue-Red 2`
#  [1] "#4A6FE3" "#5A78E2" "#6780E1" "#7489E1" "#7F91E1" "#8A99E1" "#95A2E2" "#A0AAE2" "#AAB3E2" "#B5BBE3"
#  [11] "#BFC4E3" "#C9CDE3" "#D3D5E3" "#DDDEE3" "#E3DDDE" "#E5D1D4" "#E5C6CB" "#E6BBC2" "#E6AFB9" "#E5A4B0"
#  [21] "#E498A7" "#E38D9E" "#E18195" "#DF758D" "#DC6984" "#D95C7B" "#D64E72" "#D33F6A"
#  
#  $`Blue-Red 3`
#  [1] "#002F70" "#003F84" "#19509A" "#2B60B2" "#4171C5" "#5E82CD" "#7693D5" "#8CA3DD" "#A1B3E4" "#B4C2EB"
#  [11] "#C6D0F1" "#D6DEF5" "#E5E9F8" "#F2F3F8" "#F9F2F2" "#FAE5E5" "#F7D6D6" "#F3C6C6" "#EDB4B5" "#E5A2A2"
#  [21] "#DC8F8F" "#D27B7C" "#C66767" "#B95252" "#A54142" "#8D3333" "#752425" "#5F1415"
#  
#  $`Red-Green`
#  [1] "#841859" "#9A2C6B" "#B13C7D" "#C84B90" "#D7609F" "#E376AD" "#ED8ABB" "#F59DC7" "#FDAED3" "#FFBFDE"
#  [11] "#FFCEE7" "#FFDCEF" "#FFE8F4" "#FEF2F7" "#EFF7EF" "#E0F5E0" "#D0EFD0" "#BEE8BF" "#ABDFAC" "#97D498"
#  [21] "#82C883" "#6BBB6C" "#51AD52" "#319E33" "#018E08" "#007B00" "#006900" "#005600"
#  
#  $`Purple-Green`
#  [1] "#492050" "#5D2F65" "#733E7C" "#884D93" "#9E5BAA" "#B16CBE" "#BC81C8" "#C794D1" "#D1A7DA" "#DAB8E2"
#  [11] "#E3C8E9" "#E9D6EE" "#EFE3F1" "#F1EDF2" "#EBF0EB" "#DDEADD" "#CCE3CC" "#B9D9BA" "#A5CEA5" "#8FC190"
#  [21] "#78B479" "#5FA660" "#429743" "#328533" "#287229" "#1D5F1E" "#114C12" "#023903"
#  
#  $`Purple-Brown`
#  [1] "#312A56" "#40396A" "#504880" "#5F5796" "#7066AD" "#8075C4" "#9084D9" "#A096E2" "#AFA7EA" "#BFB8F1"
#  [11] "#CEC8F7" "#DCD8FC" "#E9E7FF" "#F5F4FE" "#FDF3EE" "#FCE5D9" "#F7D5C2" "#F0C5AB" "#E6B494" "#DCA37B"
#  [21] "#D09261" "#C38044" "#AF7137" "#9A622E" "#855424" "#71451A" "#5C370E" "#492900"
#  
#  $`Green-Brown`
#  [1] "#004B40" "#005D52" "#007063" "#008274" "#009586" "#00A596" "#00B3A5" "#3EC1B4" "#69CDC2" "#89D9CF"
#  [11] "#A5E3DB" "#BEECE6" "#D6F3EF" "#EBF7F5" "#F9F3ED" "#F9EBDC" "#F4E0CA" "#EDD4B6" "#E4C6A1" "#DAB78A"
#  [21] "#CEA872" "#C19857" "#B48737" "#A2771B" "#8E670B" "#7A5700" "#664700" "#533600"
#  
#  $`Blue-Yellow 2`
#  [1] "#4F53B7" "#5E61B9" "#6C6EBD" "#797CC2" "#8788C7" "#9495CC" "#A0A1D1" "#ADAED6" "#B9BADB" "#C5C5E0"
#  [11] "#D0D0E4" "#DBDBE8" "#E5E5EC" "#EDEDEF" "#EEEEEC" "#E7E6DE" "#DEDCCF" "#D5D3BE" "#CCC8AD" "#C2BD9A"
#  [21] "#B8B286" "#ADA770" "#A39C58" "#98903C" "#8D8407" "#827900" "#776D00" "#6B6100"
#  
#  $`Blue-Yellow 3`
#  [1] "#9FA2FF" "#AAACFF" "#B3B6FF" "#BCBFFF" "#C5C7FF" "#CCCEFF" "#D3D5FF" "#D9DBFF" "#DEE0FF" "#E3E4FF"
#  [11] "#E7E8FF" "#EAEBFF" "#ECEDFF" "#EFEFFF" "#F4F1DE" "#F6F1CE" "#F6F0C2" "#F5EEB7" "#F2EBAC" "#EFE7A0"
#  [21] "#ECE394" "#E7DD87" "#E1D779" "#DBD16A" "#D4C958" "#CCC144" "#C3B827" "#BAAE00"
#  
#  $`Green-Orange`
#  [1] "#11C638" "#38C84C" "#4ECB5D" "#60CD6B" "#70D078" "#7ED285" "#8BD491" "#98D69D" "#A4D8A8" "#B0DAB3"
#  [11] "#BBDCBE" "#C7DEC8" "#D2E0D3" "#DDE2DD" "#E4DFDD" "#E7DAD2" "#E9D4C7" "#EBCEBB" "#EDC9B0" "#EFC3A4"
#  [21] "#F0BE98" "#F0B88B" "#F1B27E" "#F1AD6F" "#F1A860" "#F1A24D" "#F09D36" "#EF9708"
#  
#  $`Cyan-Magenta`
#  [1] "#0FCFC0" "#4DD2C5" "#6BD6CA" "#81D9CF" "#94DDD3" "#A4E0D8" "#B3E3DC" "#BFE5E0" "#CBE8E4" "#D5EAE7"
#  [11] "#DEECEA" "#E5EEED" "#EBEFEF" "#F0F0F0" "#F1F0F0" "#F1EEF0" "#F2EAEE" "#F2E6ED" "#F3E1EB" "#F4DBE9"
#  [21] "#F4D5E7" "#F5CEE5" "#F6C7E2" "#F6C0E0" "#F6B8DD" "#F7AFDA" "#F7A6D7" "#F79CD4"
#  
#  $Tropic
#  [1] "#009B9F" "#00A1A5" "#00A7AB" "#00ADB1" "#2EB4B6" "#53BABC" "#6CC0C2" "#82C6C8" "#95CDCF" "#A7D3D5"
#  [11] "#B8DADB" "#C9E0E1" "#D9E7E7" "#E9EDED" "#EFEBEE" "#EDE1E8" "#EAD6E3" "#E7CCDE" "#E4C1D8" "#E1B6D3"
#  [21] "#DEACCE" "#DBA1C8" "#D896C3" "#D58CBE" "#D180B9" "#CE75B4" "#CB69AF" "#C75DAA"
#  
#  $Broc
#  [1] "#002B4B" "#003658" "#004366" "#005076" "#1E5D84" "#3D6A8E" "#547898" "#6987A3" "#7E96AF" "#92A6BA"
#  [11] "#A7B6C7" "#BCC7D4" "#D3D9E1" "#EBEDF0" "#EDEDE9" "#D9D9CE" "#C7C7B6" "#B6B69E" "#A6A587" "#969571"
#  [21] "#87865A" "#787742" "#6A6926" "#5C5B00" "#4F4E00" "#424100" "#353400" "#292800"
#  
#  $Cork
#  [1] "#00294D" "#00375D" "#04456E" "#195481" "#276495" "#3D73A5" "#5882AF" "#7091B9" "#85A1C4" "#9AB0CE"
#  [11] "#AFBFD7" "#C2CEE0" "#D6DDE9" "#E8EBEF" "#E7ECE7" "#D3E0D2" "#BFD3BD" "#AAC5A8" "#95B792" "#80A97C"
#  [21] "#6A9A65" "#548C4E" "#3D7D33" "#2B6D1F" "#1F5D12" "#124D02" "#023D00" "#002F00"
#  
#  $Vik
#  [1] "#002E60" "#003B6C" "#00497C" "#00588E" "#00679F" "#1676A9" "#4485B3" "#6193BE" "#7AA3C8" "#91B2D2"
#  [11] "#A7C1DC" "#BDCFE4" "#D2DEEC" "#E6EBF1" "#EFEAE4" "#E7DACD" "#DDCAB5" "#D2BA9E" "#C6AA86" "#B9996E"
#  [21] "#AC8954" "#9F7836" "#916800" "#825800" "#714A00" "#5F3B00" "#4E2D00" "#3E2000"
#  
#  $Berlin
#  [1] "#7FBFF5" "#57ABE8" "#0D99DC" "#0087C7" "#0075AE" "#006597" "#005682" "#00486F" "#003C5D" "#00314C"
#  [11] "#00273E" "#021F30" "#091824" "#0F1316" "#171110" "#241211" "#301513" "#3D1A17" "#4B201D" "#5A2825"
#  [21] "#6C312D" "#7F3C37" "#944742" "#AB544E" "#C3625B" "#DA736C" "#E98A84" "#F8A29E"
#  
#  $Lisbon
#  [1] "#E2FCFF" "#CBEAFF" "#B3D7FB" "#9CC5ED" "#8AB3D9" "#7CA1C3" "#6F8FAD" "#627D98" "#556C83" "#485C6E"
#  [11] "#3C4B5B" "#313C47" "#262D34" "#1C1E22" "#1E1E1A" "#2D2D22" "#3C3C2B" "#4B4B34" "#5C5B3E" "#6C6C49"
#  [21] "#7D7D54" "#8F8E5F" "#A1A06A" "#B3B276" "#C5C487" "#D7D6A1" "#EAE9BA" "#FCFCD3"
#  
#  $Tofino
#  [1] "#D6E0FF" "#C3CEFF" "#B1BCF9" "#9FABE8" "#8F9AD3" "#808ABE" "#717AA9" "#626B95" "#545B81" "#464D6D"
#  [11] "#393E5A" "#2C3147" "#202435" "#161721" "#121A0F" "#1A2814" "#23371C" "#2E4625" "#39562F" "#456639"
#  [21] "#527644" "#5E8750" "#6C985C" "#79AA68" "#87BC74" "#98CD86" "#ADDE9D" "#C2EFB4"
#  
#  $ArmyRose
#  [1] "#D66982" "#D87389" "#DA7C90" "#DC8698" "#DE8F9F" "#E098A6" "#E3A2AE" "#E6ABB6" "#EAB5BE" "#EEBEC7"
#  [11] "#F3C8D0" "#F8D2D9" "#FFDDE3" "#FFE9EF" "#F3F7DA" "#E9EECB" "#DFE5BC" "#D5DCAD" "#CBD39F" "#C1C990"
#  [21] "#B8C082" "#AEB773" "#A5AF64" "#9CA655" "#929D43" "#8A943E" "#818B38" "#798233"
#  
#  $Earth
#  [1] "#A36B2B" "#A87436" "#AD7E41" "#B2874C" "#B89156" "#BD9A61" "#C2A36C" "#C7AD77" "#CDB682" "#D2C08D"
#  [11] "#D8C999" "#DED2A5" "#E4DCB0" "#EAE5BC" "#E5E6C3" "#D8DFC3" "#CCD7C2" "#C3CEC0" "#B7C7B8" "#A8C1AF"
#  [21] "#98BAA8" "#87B4A3" "#75AEA0" "#63A79E" "#509F9E" "#3F989E" "#2F8F9F" "#2686A0"
#  
#  $Fall
#  [1] "#3C5941" "#48644A" "#556E52" "#61795B" "#6E8464" "#7C8F6D" "#899B76" "#97A67F" "#A6B289" "#B4BD93"
#  [11] "#C3C99D" "#D3D4A8" "#E2E0B3" "#F2ECBE" "#F7ECBF" "#F2E2B3" "#ECD7A7" "#E8CC9A" "#E3C28D" "#DFB780"
#  [21] "#DCAB73" "#D9A066" "#D69459" "#D3884D" "#D07C42" "#CD6F39" "#CA6131" "#C7522B"
#  
#  $Geyser
#  [1] "#008585" "#2B8D86" "#449589" "#599D8B" "#6CA58F" "#7EAC94" "#8FB49B" "#9FBBA2" "#ADC3AA" "#B8CDAE"
#  [11] "#C5D5B2" "#D3DEB6" "#E2E6BB" "#F2EEC1" "#F8ECBE" "#F2E2B0" "#EED7A2" "#EACC94" "#E6C186" "#E2B678"
#  [21] "#DFAB6B" "#DC9F5E" "#D89352" "#D58747" "#D27B3D" "#CE6E35" "#CB612F" "#C7522B"
#  
#  $TealRose
#  [1] "#009593" "#0C9B91" "#3FA191" "#5BA792" "#70AC95" "#83B299" "#94B79F" "#A3BDA6" "#B1C3AE" "#BCCAB4"
#  [11] "#C5D1B7" "#D0D8BB" "#DCDFBF" "#E9E5C5" "#EDE4C5" "#E9DBBD" "#E5D3B4" "#E1CAAB" "#DEC0A3" "#DCB79A"
#  [21] "#DAAD92" "#D8A38A" "#D79984" "#D68E7F" "#D5837B" "#D47779" "#D36A78" "#D35C79"
#  
#  $Temps
#  [1] "#089392" "#019B92" "#0DA391" "#22AA8F" "#36B18D" "#4AB889" "#5DBE85" "#72C482" "#89C982" "#9DCD84"
#  [11] "#AFD287" "#C1D78B" "#D2DB91" "#E2E098" "#EADD96" "#EAD48B" "#EACA81" "#EAC079" "#EAB672" "#E9AC6D"
#  [21] "#E9A16A" "#E89768" "#E78C69" "#E5806B" "#E3756F" "#DF6A73" "#D76179" "#CF597E"
#  
#  $PuOr
#  [1] "#743700" "#894600" "#9E5500" "#B36400" "#C87300" "#DB8200" "#E99236" "#F5A259" "#FFB275" "#FFC18E"
#  [11] "#FFCFA5" "#FFDCBB" "#FFE8D1" "#FFF3E7" "#F6F6F7" "#EBEAEF" "#DDDCE6" "#CECCDD" "#BEBBD4" "#ADA9CA"
#  [21] "#9B96C0" "#8983B6" "#776FAD" "#645AA4" "#524596" "#3F347B" "#2D2262" "#1C0D51"
#  
#  $RdBu
#  [1] "#611300" "#7C1C00" "#972500" "#B32F00" "#C64200" "#D15A3C" "#DA715C" "#E18778" "#E79C91" "#EBB0A8"
#  [11] "#EFC4BE" "#F2D6D2" "#F4E6E5" "#F7F4F4" "#F7F8F8" "#ECF1F5" "#DDE8F0" "#C9DDEA" "#B3D0E4" "#98C2DC"
#  [21] "#79B3D4" "#4FA4CB" "#0093C1" "#0083B8" "#0073B0" "#005E94" "#004978" "#003560"
#  
#  $RdGy
#  [1] "#68001D" "#810D21" "#9A1D20" "#B32B16" "#CC3B00" "#D55827" "#DD7145" "#E38860" "#E99D7A" "#EDB193"
#  [11] "#F1C5AC" "#F3D6C5" "#F6E7DC" "#F8F4F1" "#F5F5F5" "#E9E9E9" "#DCDCDC" "#CDCDCD" "#BEBEBE" "#AFAFAF"
#  [21] "#9F9F9F" "#8F8F8F" "#7F7F7F" "#6F6F6F" "#5F5F5F" "#4F4F4F" "#404040" "#303030"
#  
#  $PiYG
#  [1] "#90005D" "#AA0070" "#C40983" "#CF3C91" "#D75A9E" "#DF73AA" "#E689B7" "#EB9EC2" "#F0B1CE" "#F4C2D8"
#  [11] "#F7D2E2" "#F9E1EB" "#FAEDF2" "#FAF7F8" "#F3FAF0" "#E7F8DE" "#DAF2CC" "#CBEBB9" "#BBE2A4" "#AAD88D"
#  [21] "#99CD75" "#86C059" "#73B336" "#63A312" "#579209" "#4B8100" "#3E6F00" "#315D00"
#  
#  $PRGn
#  [1] "#410E48" "#551D5D" "#692B73" "#7F388B" "#92489F" "#9D60A8" "#A876B1" "#B38ABB" "#BE9EC5" "#C9B1CF"
#  [11] "#D5C3D9" "#DFD4E2" "#EAE4EB" "#F3F1F3" "#F0F7F0" "#E2F3E3" "#D2EDD3" "#C0E5C1" "#ADDBAD" "#98CF99"
#  [21] "#82C282" "#69B46A" "#51A552" "#479248" "#3C7F3D" "#316C31" "#245825" "#174417"
#  
#  $BrBG
#  [1] "#533600" "#674600" "#7B5700" "#906700" "#A47700" "#B98600" "#C7973E" "#D4A760" "#DFB67C" "#E9C596"
#  [11] "#F1D3AD" "#F7DFC4" "#FBEAD8" "#FAF3EC" "#EEF6F5" "#DDF0ED" "#CAE8E3" "#B6DDD8" "#A0D1CB" "#88C4BC"
#  [21] "#6EB6AD" "#50A79D" "#26978D" "#00877C" "#00746A" "#006158" "#004F46" "#003D34"
#  
#  $RdYlBu
#  [1] "#A51122" "#B92400" "#CC3600" "#DC4A00" "#E16400" "#E67A00" "#E98E00" "#ECA105" "#EFB343" "#F1C363"
#  [11] "#F4D27E" "#F7E195" "#FAEDA9" "#FDF9B9" "#F4FACE" "#E5F0D6" "#DBE4D7" "#C6DCC7" "#ACD2BB" "#90C9B2"
#  [21] "#71BEAC" "#4EB2A9" "#20A5A8" "#0097A8" "#0088A7" "#0077A6" "#0063A4" "#324DA0"
#  
#  $RdYlGn
#  [1] "#A51122" "#B92400" "#CC3600" "#DC4A00" "#E16400" "#E67A00" "#E98E00" "#ECA105" "#EFB343" "#F1C363"
#  [11] "#F4D27E" "#F7E195" "#FAEDA9" "#FDF9B9" "#F9FAB5" "#EAF3A3" "#DAEA91" "#C8E07E" "#B4D66C" "#9FCB59"
#  [21] "#8AC048" "#76B345" "#63A642" "#4F983E" "#3A8B3A" "#207D35" "#00702F" "#006228"
#  
#  $Spectral
#  [1] "#A71B4B" "#B82F47" "#C9413F" "#D8522F" "#E66407" "#EB7A04" "#F08D17" "#F4A02E" "#F7B245" "#F9C25C"
#  [11] "#FBD273" "#FCE08A" "#FEEEA0" "#FEF9B5" "#F4FCBA" "#DBF6B3" "#BFEFAF" "#A1E7AD" "#81DEAD" "#5DD3B0"
#  [21] "#31C7B2" "#00B9B4" "#00ABB6" "#009BB6" "#008AB4" "#1E77AF" "#4262A8" "#584B9F"
#  
#  $`Zissou 1`
#  [1] "#3B99B1" "#35A0AE" "#37A6A9" "#40ABA5" "#4EB0A0" "#5EB49C" "#6FB798" "#80BB95" "#90BD94" "#9FC095"
#  [11] "#ABC391" "#B8C68E" "#C7C98A" "#D8CB81" "#EAC628" "#EABE23" "#EAB51F" "#E9AC1C" "#E8A419" "#E79B14"
#  [21] "#E7910D" "#E78802" "#E77D00" "#E97100" "#EA6400" "#ED5400" "#F03F08" "#F5191C"
#  
#  $Cividis
#  [1] "#00214E" "#002757" "#002D61" "#00336B" "#18396D" "#29406D" "#36476D" "#414E6E" "#4A5570" "#545C71"
#  [11] "#5D6373" "#666A75" "#6F7178" "#78787A" "#81807C" "#8B877B" "#958F79" "#9F9778" "#A89F76" "#B2A773"
#  [21] "#BCAF70" "#C6B76D" "#D0BF69" "#DAC863" "#E4D05D" "#EED855" "#F9E14C" "#FFE93F"
#  
#  $Roma
#  [1] "#7D0112" "#8B2B0E" "#99440A" "#A55A0B" "#B06F15" "#B98223" "#C09533" "#C6A745" "#CAB858" "#CCC86C"
#  [11] "#CDD680" "#CDE294" "#CDEBA7" "#CEF1BC" "#CEEDCB" "#B9E5C3" "#A1DCBB" "#86D2B5" "#67C7B1" "#42BBAE"
#  [21] "#00AEAB" "#00A0A9" "#0091A7" "#0081A5" "#0070A3" "#005DA0" "#00489F" "#1F28A2"
#-----myFunctions------------
mypreprocess = function(count,meta.data=NULL,min.cells=3,min.features=200,res=0.5,npcs=50,dim=15){
    #-----------S1-data preprocess -------------------------------------------
    if(!is.null(meta.data=meta.data)){
       pbmc <- CreateSeuratObject(count,meta.data =meta.data ,min.cells = min.cells, min.features = min.features)
    }else{
       pbmc <- CreateSeuratObject(count,min.cells = min.cells, min.features = min.features)
       # pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
       pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)
       # pbmc <- NormalizeData(pbmc)
    }
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    pbmc <-SCTransform(pbmc)
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    pbmc <- FindVariableFeatures(pbmc, x.low.cutoff = 0.0125, y.cutoff = 0.25, do.plot=FALSE)
    pbmc <- RunPCA(pbmc,features = c(VariableFeatures(pbmc),'mCherry-CAR'),npcs = npcs)
    pbmc <- FindNeighbors(pbmc,dim=1:dim)
    pbmc <- FindClusters(pbmc,dim=1:dim,resolution = res)
    pbmc <- RunUMAP(pbmc,dim=1:dim)
    pbmc <- RunTSNE(pbmc)
 }

mm.genes= readxl::read_xlsx("~/LJX/mm/Cell_marker_Mouse.xlsx")
ref_celltype = mm.genes[,c("Symbol","cell_name")]

myAUCell = function(exprMatrix,geneSets,plotStats=FALSE,plot=F,auc.threshold=0.3){
   cells_AUC <- AUCell_run(exprMatrix, geneSets)
   cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=plotStats)
   cells_AUC <- AUCell_run(exprMatrix, geneSets, aucMaxRank=nrow(cells_rankings)*0.05,)
   
   if(plot==TRUE){
      set.seed(12345)
      par(mfrow=c(2,3))
      cells_assignment <- AUCell_exploreThresholds(cells_AUC,plotHist = TRUE,assign=T )
      selectedThresholds <- getThresholdSelected(cells_assignment)
     # Splits the plot into two rows and three columns
      for(geneSetName in names(selectedThresholds)[selectedThresholds<auc.threshold])
      {
         nBreaks <- 5 # Number of levels in the color palettes
         # Color palette for the cells that do not pass the threshold
         colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
         # Color palette for the cells that pass the threshold
         colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
         
         # Split cells according to their AUC value for the gene set
         passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
         if(sum(passThreshold) >0 )
         {
            aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
            
            # Assign cell color
            cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                           setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
            
            # Plot
            plot(cellsTsne, main=geneSetName,
                 sub="Pink/red cells pass the threshold",
                 col=cellColor[rownames(cellsTsne)], pch=16) 
         }
      }
      
      return(cells_AUC)
      
   }else{
      
      return(cells_AUC)
   }
   
} 
myenrichrCelltype =function(obj,ref_celltype,tc.mega,logfc=0,n=5){
   markers = FindAllMarkers(obj,only.pos = T,logfc.threshold = logfc)
   markers = markers[markers$pct.1>n*markers$pct.2,]
   genelist = split(markers$gene,markers$cluster)
   celltype = lapply(1:length(genelist),function(i){
      celltype.enrichr = enricher(gene = genelist[[i]],
                                  TERM2GENE = ref_celltype,pvalueCutoff = 1)
      if(!is.null(celltype.enrichr)){
         celltype.enrichr = celltype.enrichr@result
         # celltype.enrichr$cluster = rep(names(genelist)[i],nrow(celltype.enrichr))
         return(celltype.enrichr)
      }
      
   })
   
   names(celltype) = names(genelist)
   
   celltype1 =data.table::rbindlist(celltype,use.names = T,fill = T,idcol = "subcluster")
   
   # celltype1$MegaTC = rep(tc.mega,nrow(celltype1))
   return(celltype1)
}#enrichr for Cell Marker annotation
myenrichr =function(markers,ref_celltype = ref_celltype){
   genelist = split(markers$gene,markers$cluster)
   celltype = lapply(1:length(genelist),function(i){
      celltype.enrichr = enricher(gene = genelist[[i]],TERM2GENE = ref_celltype,pvalueCutoff = 1)
      if(!is.null(celltype.enrichr)){
         celltype.enrichr = celltype.enrichr@result
         # celltype.enrichr$cluster = rep(names(genelist)[i],nrow(celltype.enrichr))
         return(celltype.enrichr)
      }
      
   })
   
   names(celltype) = names(genelist)
   
   celltype1 =data.table::rbindlist(celltype,use.names = T,fill = T,idcol = "subcluster")
}
#st cluster/pos/distance
cal_zoom_rate = function(width, height){
   std_width = 1000
   std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
   if(std_width / std_height > width / height){
      scale = width / std_width
   }
   else{
      scale = height / std_height
   }
   return(scale)
}

#-----create st seurat obj ---------
mySeurat = function(path,levelpath,hefigpath,
                    projname){
   setwd(path)
   cd = Read10X(levelpath)
   setwd(levelpath)
   pos = readr::read_delim(gzfile("barcodes_pos.tsv.gz"),col_names = F)
   pos= data.frame(pos,row.names = 1)
   colnames(pos)=c("x","y")
   
   he_fig =png::readPNG(hefigpath)
   
   zoom_scale = cal_zoom_rate(ncol(he_fig), nrow(he_fig))
   pos$x=pos$x * zoom_scale
   pos$y=pos$y * zoom_scale
   pos$y=dim(he_fig)[1]-pos$y
   
   # setwd("..")
   data = CreateSeuratObject(cd,project = projname,min.cells = 3)
   data= PercentageFeatureSet(data,pattern = "^mt-",col.name = "percent.mt")
   data = SCTransform(data,return.only.var.genes = F,
                      # method="glmGamPoi",
                      vars.to.regress = "percent.mt",min_cells=3)
   data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
   # data = subset(data, subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mt < 5)
   all.genes <- rownames(data)
   data <- ScaleData(data, features = all.genes)
   # data = RunPCA(data)
   # npc = 20
   # data = basicFindCluster(data,npc=npc)
   feature =  intersect(c(VariableFeatures(data),features,"EGFP","mCherry-CAR"),rownames(data@assays$RNA@scale.data))%>%unique()
   data = AddMetaData(data,pos)
   data = RunPCA(data,
                 features =feature )
   data =FindNeighbors(data,dim=1:30)
   data = FindClusters(data,resolution = 1)
   data = RunUMAP(data,dim=1:30)
   data = RunTSNE(data)
   return(data)
}


myclstplot = function(he_file,FilePath,
                      Cluster,hetype="png",
                      celltypecol = NULL,pt.size=3.5){
   
   if(hetype=="png"){
      he_fig=png::readPNG(he_file)
   }if(hetype=='tiff'){
      he_fig=tiff::readTIFF(he_file)
   }
   
   w=ncol(he_fig)
   h=nrow(he_fig)
   
   cal_zoom_rate = function(width, height){
      std_width = 1000
      std_height = std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
      if(std_width / std_height > width / height){
         scale = width / std_width
      }
      else{
         scale = height / std_height
      }
      return(scale)
   }
   
   zoom_scale = cal_zoom_rate(ncol(he_fig), nrow(he_fig))
   bc_pos_file = gzfile(paste(FilePath,"barcodes_pos.tsv.gz", sep  = "/"),'rt')
   bc_pos = read.table(bc_pos_file, header = FALSE, sep = '\t', quote = '')
   names(bc_pos)=c("Barcode","x","y")
   Cluster = data.frame(cluster=Cluster)
   colnames(Cluster)=c("cluster")
   Cluster$Barcode = rownames(Cluster)
   bc_pos=merge(bc_pos,Cluster,by="Barcode",all=FALSE)
   
   bc_pos$x=bc_pos$x * zoom_scale
   bc_pos$y=bc_pos$y * zoom_scale
   bc_pos$y=dim(he_fig)[1]-bc_pos$y
   
   clstr_plot<-ggplot(data=bc_pos,aes(x ,y)) + 
      ggpubr::background_image(he_fig)+
      # geom_point(aes(colour = cluster),size = opt$point_size, shape=16) + 
      geom_point(aes(fill = cluster,color=cluster),size = pt.size, shape=21) + 
      # scale_fill_brewer(palette = "Set1")+
      # scale_color_manual(values = col)+
      theme_bw() +
      theme(plot.title = element_text(face = 2,size = 50,hjust = 0.5)) +
      theme(axis.ticks = element_blank(), 
            axis.text.y = element_blank(),
            panel.border = element_blank(), 
            axis.text.x = element_blank()) +
      xlab('')+ylab('')+
      coord_cartesian(xlim = c(0, dim(he_fig)[2]), ylim = c(0, dim(he_fig)[1]), 
                      expand = FALSE)+
      guides(colour = guide_legend(override.aes = list(size=3.5)))+
      theme(legend.position = "none")
   
   # if(is.null(celltypecol)){
   #   clstr_plot = clstr_plot
   # }else{
   #   clstr_plot = clstr_plot+scale_fill_manual(values= celltypecol)+
   #     scale_color_manual(values= celltypecol)
   # }
   # 
   return(clstr_plot)
} 

st_distance = function(x,y){
   central_x = (max(x)-min(x))/2
   central_y =(max(y)-min(y))/2
   distance = sqrt((x-central_x)^2 + (y-central_y)^2)
}

myAddspatial = function(pbmc){
  pbmc@reductions$spatial <- pbmc@reductions$umap
  spatial_embedding = data.frame(pbmc@meta.data[,c("x","y")])
  colnames(spatial_embedding)=c("SPATIAL_1","SPATIAL_2")
  pbmc@reductions$spatial@cell.embeddings <- as.matrix(spatial_embedding)
  pbmc@reductions$spatial@key="SPATIAL_"
  return(pbmc)
}


cd8tmarkers

markers =c('mCherry-CAR','Cd3e',
           # 'Cd2','Cd5','Cd7',
           'Cd8b1',
           'Ccr7','Sell',
           'Itgae','Nkg7',
           'Tcf7','Cd44',
           'Il2ra','Il2rb',
           'Tgfb1','Il10',
           'Cd28','Lck','Tnfrsf9',
           'Gzma','Gzmb','Prf1','Tnf','Ifng',
           'Tox','Pdcd1','Tight','Ctla4','Pseudotime')

mydot(cd8t,markers)|DimPlot(cd8t,label = T,repel = T)
mydot(cd8t,markers)|FeaturePlot(cd8t,features = 'Pseudotime',label = T,repel = T)
distance.list = lapply(1:length(sc_all_list),function(i){
  x=sc_all_list[[i]]@spatial_location$x
  y=sc_all_list[[i]]@spatial_location$y
  distance = st_distance(x,y)
  return(distance)
})


df = sc_all_list[[1]]@Proportion_CARD
df =reshape2::melt(df)
df$distance = rep(distance.list[[1]],
                  levels(df$Var2) %>%length())


p1=FeaturePlot(st.list$Nonresponder,features = c("mCherry-CAR","Tox","Ctla4","Pdcd1"),ncol=2,
               reduction = "spatial",
               cols = rcolors::rcolors$MPL_Reds) & xlim(0,1000) &ylim(0,1000)
p2=FeaturePlot(st.list$Responder,features = c("mCherry-CAR","Tox","Ctla4","Pdcd1"),ncol=2,
               reduction = "spatial",
               cols = rcolors::rcolors$MPL_Reds) & xlim(0,1000) &ylim(0,1000)

#======Fig panel CARD deconvolution ===========================================
 library(Seurat)
 library(CARD) #NBT,2022
 library(dplyr)
 
 n9_path = "~/LJX/st/n9/05.AllheStat/heAuto_level_matrix/subdata"
 r5_path ="~/LJX/st/r5/05.AllheStat/heAuto_level_matrix/subdata"
 levelpath = dir(n9_path)
 n9_fig_path = '~/LJX/st/n9/05.AllheStat/allhe/he_roi_small.png'
 r5_fig_path = '~/LJX/st/r5/05.AllheStat/allhe/he_roi_small.png'
 # batchcolor = c("#1f77b4", "#ff7f0e")
 # pycolor2=c("#1f77b4","#ff7f0e" ,"#279e68","#d62728","#aa40fc" ,"#8c564b" ,"#e377c2","#b5bd61" ,"#17becf","#aec7e8","#ffbb78" ,"#98df8a") 
 # cellname=c("Fibroblast" , "DNT","NKT" , "CD4T Cell" , "Myeloid","CD8T Cell", "DPT" ,"T Cell MKI67+" ,"CAFs","Tumor","B Cell", "Epithelial cell" )
 # names(pycolor2)=cellname
 celltypecolor =rcolors::rcolors$srip_reanalysis[2:19]
 names(celltypecolor)=levels(immune.combined@active.ident)
 st.n9 = mySeurat()
 
 
#======scRNA referenced celltype deconvolution================================
#--------DATA Preparation---------------
# 
#  load("./spatial_count.RData")
#  spatial_count[1:4,1:4]
#  

#  
#  load("./spatial_location.RData")
#  spatial_location[1:4,]
#  x  y
#  10x10 10 10
#  10x13 10 13
#  10x14 10 14
#  10x15 10 15
#  load("./sc_meta.RData")
#  
#  sc_meta[1:4,]
#  cellID                             cellType sampleInfo
#  Cell1  Cell1                         Acinar_cells    sample1
#  Cell2  Cell2          Ductal_terminal_ductal_like    sample1
#  Cell3  Cell3          Ductal_terminal_ductal_like    sample1
#  Cell4  Cell4 Ductal_CRISP3_high-centroacinar_like    sample1

 
# setwd("~/LJX/st/L13_cluster/CARD_res/")
# load("~/LJX/sc_in/adata_preprocess.RData")
# # sc_meta = readr::read_csv("~/LJX/sc_in/adata_sc_preprocessed.csv")

# adata
# cellsplit = split(names(adata@active.ident),adata$sample)
# sc.list = lapply(cellsplit,function(x){
#    obj = subset(adata,cells=x)
# })
# 
# sc.list = list("sc_all" = adata,"myesub"=mye)

# rm(adata)

setwd("~/LJX/sc_in/")
n9 = Read10X("./N9/")
r5 = Read10X("./R5/")
obj.list = list('Nonresponder'=n9,'Responder'=r5)
obs = readr::read_csv("obs/obs_all_20230815.csv")

n9.cd = obj.list$`Non-response`[,obs$barcodeID[obs$sample=="N9"]]
r5.cd = obj.list$Response[,obs$barcodeID[obs$sample=="R5"]]

colnames(n9.cd)=paste0(colnames(n9.cd),"-1")
colnames(r5.cd)=paste0(colnames(r5.cd),"-2")

barcode =c(colnames(n9.cd),colnames(r5.cd))
barcode2 = ifelse(obs$sample=="N9",paste0(obs$barcodeID,"-1"),paste0(obs$barcodeID,"-2"))

all.equal(barcode,barcode2)

sc_count = cbind(n9.cd,r5.cd)
sc_count =sc_count[,barcode2]
is.na(colnames(sc_count)) %>% table()

# sc_meta = data.frame(cellID=barcode2,
#                      
#                      cellType=obs$Cell_type,
#                      sampleInfo = obs$sample)
# sc_meta$cellType = factor(sc_meta$cellType,levels = names(pycolor2))

# sc_meta = data.frame(cellID=names(adata@active.ident),cellType=adata$Cell_type,sampleInfo = adata$sample)
# sc_count =adata@assays$RNA@counts


rownames(sc_meta)=sc_meta$cellID
pathlist = list(n9_path,r5_path)

colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")

ctcolors =c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
            '#7f7f7f', '#bcbd22', '#17becf',"#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
            "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
            "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")

# dir.create("Deconvoluted_to_sc_all")
# setwd("Deconvoluted_to_sc_all")
ctcolors = pycolor2
ctcolors=celltypecolor
names(pathlist)=c("Nonresponder","Responder")
sc_count = immune.combined@assays$RNA@counts
# sc_meta = immune.combined@meta.data
sc_meta = data.frame(cellID=names(immune.combined@active.ident),
                     # cellType=immune.combined@active.ident,
                     cellType=sc_meta$cellType,
                     cluster = immune.combined$seurat_clusters,
                     sampleInfo = immune.combined$sample)
st.list = SplitObject(st.integrated,split.by = 'sample')

sc_all_list = lapply(1:2,function(i){
   path0 = pathlist[[i]]
   setwd(path0)
   # levelpath = dir(path = path0)
   # levelpath = grep("pdf",levelpath,invert = T,value = T)
   
   spatial_count =  st.list[[i]]@assays$RNA@counts
   # setwd(levelpath[2])
   # spatial_location=readr::read_delim(gzfile("barcodes_pos.tsv.gz"),col_names = F)
   spatial_location = st.list[[i]]@reductions$spatial@cell.embeddings
   

   # spatial_location = data.frame(x=spatial_location$x,y=spatial_location$y)
   # rownames(spatial_location)=rownames(st.list[[i]]@meta.data)
   colnames(spatial_location)=c("x","y")
   
   setwd("~/LJX/st/L13_cluster/")
   # dir.create("CARD_res")
   setwd("./CARD_res/Deconvoluted_to_sc_all")
   sc_count =GetAssayData(reference)
   sc_meta = reference@meta.data
   # hcl.colors()
   #-------create CARD OBJ and deconvolution-----------
   CARD_obj = createCARDObject(
      sc_count = sc_count,
      sc_meta = sc_meta,
      spatial_count = spatial_count,
      spatial_location = spatial_location,
      ct.varname = "celltype",
      ct.select = unique(sc_meta$cellType),
      sample.varname = 'sample',
      minCountGene = 100,
      minCountSpot = 3) 
   
   CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
   CARD_obj = CARD.imputation(CARD_obj,NumGrids = 5000,ineibor = 10,exclude = NULL)
   return(CARD_obj)
   
   # #----------plot results---------------------------------------------------------
   # pdf(paste0("deconv_all_",names(pathlist)[i],".pdf"))
   # p1 <- CARD.visualize.pie(
   #    proportion = CARD_obj@Proportion_CARD,
   #    spatial_location = CARD_obj@spatial_location,
   #    colors = ctcolors,
   #    radius = 7)
   # 
   # print(p1)
   # 
   # ## select the cell type that we are interested
   # # ct.visualize = c("Acinar_cells","Cancer_clone_A","Cancer_clone_B","Ductal_terminal_ductal_like","Ductal_CRISP3_high-centroacinar_like","Ductal_MHC_Class_II","Ductal_APOL1_high-hypoxic","Fibroblasts")
   # ct = levels(sc_meta$cellType)
   # ct.visualize=ct
   # ## visualize the spatial distribution of the cell type proportion
   # 
   # p2 <- CARD.visualize.prop(
   #    proportion = CARD_obj@Proportion_CARD,        
   #    spatial_location = CARD_obj@spatial_location, 
   #    ct.visualize = ct.visualize,                 ### selected cell types to visualize
   #    colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
   #    NumCols = 4,                                 ### number of columns in the figure panel
   #    pointSize = 0.1)                             ### point size in ggplot2 scatterplot  
   # 
   # print(p2)
   # 
   # ## visualize the spatial distribution of two cell types on the same plot
   # p3 = CARD.visualize.prop.2CT(
   #    proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
   #    spatial_location = CARD_obj@spatial_location,                      ### spatial location information
   #    ct2.visualize = c("CD4T Cell"  ,   "CD8T Cell"  ),              ### two cell types you want to visualize
   #    colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))       ### two color scales                             
   # 
   # print(p3)
   # 
   # p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
   # 
   # print(p4)
   # 
   # CARD_obj = CARD.imputation(CARD_obj,NumGrids = 5000,ineibor = 10,exclude = NULL)
   # 
   # location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
   #                                        y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
   # rownames(location_imputation) = rownames(CARD_obj@refined_prop)
   # library(ggplot2)
   # p5 <- ggplot(location_imputation, 
   #              aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
   #    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
   #          legend.position="bottom",
   #          panel.background = element_blank(),
   #          plot.background = element_blank(),
   #          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
   # print(p5)
   # 
   # 
   # p6 <- CARD.visualize.prop(
   #    proportion = CARD_obj@refined_prop,                         
   #    spatial_location = location_imputation,            
   #    ct.visualize = ct.visualize,   pointSize = 0.5,                 
   #    colors = c("lightblue","lightyellow","red"),    
   #    NumCols = 4)                                  
   # print(p6)
   # 
   # features = c('Lyz2',"Cd14","Adgre1" ,"Gzmb","Fcgr3","Klrk1","Mki67","Il2ra","Cd8a","Col1a1")
   # features =intersect(features,rownames(CARD_obj@refined_expression))
   # p7 <- CARD.visualize.gene(
   #    spatial_expression = CARD_obj@refined_expression,
   #    spatial_location = location_imputation,
   #    gene.visualize = features,
   #    colors = NULL,
   #    NumCols = 6)
   # print(p7)
   # 
   # features = c("Cd14","Adgre1" ,"Gzmb","Fcgr3","Klrk1","Mki67","Il2ra","Cd8a","Ly6a","Col1a1")
   # features =intersect(features,rownames(CARD_obj@spatial_countMat))
   # p8 <- CARD.visualize.gene(
   #    spatial_expression = CARD_obj@spatial_countMat,
   #    spatial_location = CARD_obj@spatial_location,
   #    gene.visualize = features,
   #    colors = NULL,
   #    NumCols = 6)
   # print(p8)
   # 
   # dev.off()
   # return(CARD_obj)
   # 
})


p1.list = lapply(1:length(sc_all_list),function(x){
  CARD_obj = sc_all_list[[x]]
  p1 <- CARD.visualize.pie(
    proportion = CARD_obj@Proportion_CARD,
    spatial_location = CARD_obj@spatial_location,
    # colors = ctcolors,
    radius = 7)
  return(p1)+ggtitle(names(sc_all_list)[x])& xlim(0,1000) & ylim(0,1000)
}) 

n9.proportion_card = sc_all_list[[1]]@Proportion_CARD
n9.pos = st.list$Nonresponder@meta.data
r5.proportion_card = sc_all_list[[2]]@Proportion_CARD
r5.pos = st.list$Responder@meta.data

df = cbind(sc_all_list[[1]]@Proportion_CARD,distance.list[[1]])

df = as.data.frame(df)
colnames(df)=c(colnames(df)[1:28],"distance")
df = df[order(df$distance),]
pheatmap::pheatmap(t(df[,!colnames(df)%in%c("distance")]),
                   show_colnames = F,cluster_cols = F,
                   annotation_col = data.frame(Distance=distance.list[[1]],
                                               Tumor = df$Tumor,
                                               row.names = rownames(sc_all_list[[1]]@Proportion_CARD)),
                   scale = "row",cutree_rows = 3,border_color = NA,
                   color = colorRampPalette(rcolors::rcolors$temp_19lev)(50),treeheight_row = 0,
                   filename ="nonresopnder.pdf" ,width = 6,height = 4)

ggsave("p1_CARD_all.png",p1.list[[1]]|p1.list[[2]],width = 12,height = 7)

markers = FindAllMarkers(immune.combined,only.pos = T)
markers =markers[markers$pct.1>markers$pct.2,]
genelist = split(markers$gene,markers$cluster)
gores =lapply(genelist,function(x){
  go=enrichGO(x,OrgDb = org.Mm.eg.db,
              keyType = "SYMBOL",pvalueCutoff = 0.05,ont = "BP")
  return(go@result)
})
gores = data.table::rbindlist(gores,fill = T,use.names = T,idcol = "cellType")

gores$cellType =factor(gores$cellType,
                       levels=unique(tc_go$cellType)[c(21,22,8,14,9,5,25,10,28,1:3,13,11,6,27,17,12,24,23,4,26,7,15,18:20,16)])
tc_go = grep("T cell",gores$Description,value = T) 
tc_go =table(tc_go)%>%as.data.frame()
tc_go = tc_go[tc_go$Freq<14,]
tc_go =gores[gores$Description%in%tc_go$tc_go,]
p =ggplot(tc_go,aes(x=cellType,y=reorder(Description,as.numeric(cellType)),
                    color=-log10(pvalue),size=Count))+
  geom_point()+theme_bw()+xlab("")+ylab("")+scale_size_continuous(range = c(3,6))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
p=p+scale_color_viridis_c()
ggsave("Function_Enrichment.pdf",p,width = 10,height =8)

ct.visualize = levels(immune.combined@active.ident)

p2.list = lapply(1:length(sc_all_list),function(x){
  CARD_obj = sc_all_list[[x]]
  p2 <- CARD.visualize.prop(
    proportion = CARD_obj@Proportion_CARD,
    spatial_location = CARD_obj@spatial_location,
    ct.visualize = ct.visualize,                 ### selected cell types to visualize
    colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
    NumCols = 4,                                 ### number of columns in the figure panel
    pointSize = 0.1)                             ### point size in ggplot2 scatterplot
  
  return(p2+ggtitle(names(sc_all_list)[x])& xlim(0,1000) & ylim(0,1000))
  
  
})

ggsave("p2_CARD_all.png",p2.list[[1]]|p2.list[[2]],width = 24,height = 24)




#------CARD BY SEURAT CLUSTERS ------------------------------------
clustercolor = hcl.colors(palette = "Set 2",
                          n = length(levels(immune.combined$seurat_clusters)))
names(clustercolor)=levels(immune.combined$seurat_clusters)
# clustercolor
# "#ED90A4" "#EB9397" "#E7968A" "#E39A7D" "#DD9D71" "#D5A166" "#CDA55C" "#C4A954" "#B9AD50"
# "#AEB050" "#A1B454" "#93B75B" "#84BA64" "#73BC6F" "#60BE7B" "#4AC087" "#2EC194" "#00C1A0"
# "#00C1AC" "#00C1B8" "#00C0C3" "#00BECD" "#1FBBD5" "#44B8DD" "#5FB4E3" "#77B0E8" "#8CABEB"
# "#9FA6EC" "#B0A1EC" "#BE9CEA" "#CB98E6" "#D694E0" "#DF91D9" "#E58FD0" "#EA8EC6" "#ED8EBB"
# "#EE8FB0"
sc_cluster_list = lapply(1:2,function(i){
   path0 = pathlist[[i]]
   setwd(path0)
   levelpath = dir(path = path0)
   levelpath = grep("pdf",levelpath,invert = T,value = T)
   
   spatial_count =  st.list[[i]]@assays$RNA@counts
   # setwd(levelpath[2])
   # spatial_location=readr::read_delim(gzfile("barcodes_pos.tsv.gz"),col_names = F)
   spatial_location = st.list[[i]]@meta.data
   spatial_location = data.frame(x=spatial_location$x,y=spatial_location$y)
   rownames(spatial_location)=rownames(st.list[[i]]@meta.data)
   colnames(spatial_location)=c("x","y")
   
   setwd("~/LJX/st/L13_cluster/")
   # dir.create("CARD_res")
   setwd("./CARD_res/Deconvoluted_to_sc_all")
   
   # hcl.colors()
   #-------create CARD OBJ and deconvolution-----------
   CARD_obj = createCARDObject(
      sc_count = sc_count,
      sc_meta = sc_meta,
      spatial_count = spatial_count,
      spatial_location = spatial_location,
      ct.varname = "cluster",
      ct.select = unique(sc_meta$cluster),
      sample.varname = "sampleInfo",
      minCountGene = 100,
      minCountSpot = 3) 
   
   CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
   return(CARD_obj)

})








#-------- CARD--mye for subset----------------------------------
adata=mye
cellmeta = readr::read_csv('~/LJX/sc_in/metafile/obs_mye_sub_npc15_res1.csv')
sc_meta = data.frame(cellID=names(mye@active.ident),cellType=cellmeta$mye_celltype,sampleInfo = cellmeta$sample)
sc_count =adata@assays$RNA@counts
rownames(sc_meta)=sc_meta$cellID
sc_meta$cellType= factor(sc_meta$cellType,levels=c("Mono","M1","TAM","cDC","pDC"))

sc_mye_list = lapply(1:2,function(i){
   path0 = pathlist[[i]]
   setwd(path0)
   levelpath = dir(path = path0)
   levelpath = grep("pdf",levelpath,invert = T,value = T)
   
   spatial_count =  Read10X(levelpath[2])
   setwd(levelpath[2])
   spatial_location=readr::read_delim(gzfile("barcodes_pos.tsv.gz"),col_names = F)
   spatial_location = data.frame(spatial_location,row.names = 1)
   colnames(spatial_location)=c("x","y")
   
   setwd("~/LJX/st/L13_cluster/")
   # dir.create("CARD_res")
   setwd("./CARD_res/Deconvoluted_to_sc_all")
   
   # hcl.colors()
   #-------create CARD OBJ and deconvolution-----------
   CARD_obj = createCARDObject(
      sc_count = sc_count,
      sc_meta = sc_meta,
      spatial_count = spatial_count,
      spatial_location = spatial_location,
      ct.varname = "cellType",
      ct.select = unique(sc_meta$cellType),
      sample.varname = "sampleInfo",
      minCountGene = 100,
      minCountSpot = 3) 
   
   CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
   
   #----------plot results---------------------------------------------------------
   pdf(paste0("deconv_mye_",names(pathlist)[i],".pdf"))
   p1 <- CARD.visualize.pie(
      proportion = CARD_obj@Proportion_CARD,
      spatial_location = CARD_obj@spatial_location, 
      colors = ctcolors, 
      radius = 5)
   
   print(p1)
   
   
   ## select the cell type that we are interested
   # ct.visualize = c("Acinar_cells","Cancer_clone_A","Cancer_clone_B","Ductal_terminal_ductal_like","Ductal_CRISP3_high-centroacinar_like","Ductal_MHC_Class_II","Ductal_APOL1_high-hypoxic","Fibroblasts")
   ct = levels(sc_meta$cellType)
   ct.visualize=ct
   ## visualize the spatial distribution of the cell type proportion
   
   p2 <- CARD.visualize.prop(
      proportion = CARD_obj@Proportion_CARD,        
      spatial_location = CARD_obj@spatial_location, 
      ct.visualize = ct.visualize,                 ### selected cell types to visualize
      colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
      NumCols = 4,                                 ### number of columns in the figure panel
      pointSize = 0.1)                             ### point size in ggplot2 scatterplot  
   
   print(p2)
   
   ## visualize the spatial distribution of two cell types on the same plot
   p3 = CARD.visualize.prop.2CT(
      proportion = CARD_obj@Proportion_CARD,                             ### Cell type proportion estimated by CARD
      spatial_location = CARD_obj@spatial_location,                      ### spatial location information
      ct2.visualize = c("M1"  ,   "TAM"  ),              ### two cell types you want to visualize
      colors = list(c("lightblue","lightyellow","red"),c("lightblue","lightyellow","black")))       ### two color scales                             
   
   print(p3)
   
   p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD,colors = NULL) # if not provide, we will use the default colors
   
   print(p4)
   
   CARD_obj = CARD.imputation(CARD_obj,NumGrids = 5000,ineibor = 10,exclude = NULL)
   
   location_imputation = cbind.data.frame(x=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",1)),
                                          y=as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop),split="x"),"[",2)))
   rownames(location_imputation) = rownames(CARD_obj@refined_prop)
   library(ggplot2)
   p5 <- ggplot(location_imputation, 
                aes(x = x, y = y)) + geom_point(shape=22,color = "#7dc7f5")+
      theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            legend.position="bottom",
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
   print(p5)
   
   
   p6 <- CARD.visualize.prop(
      proportion = CARD_obj@refined_prop,                         
      spatial_location = location_imputation,            
      ct.visualize = ct.visualize,   pointSize = 0.5,                 
      colors = c("lightblue","lightyellow","red"),    
      NumCols = 4)                                  
   print(p6)
   
   features = c("Ly6g","Itgam","Itgax","Sigel","Cd14","Adgre1" ,"Cd163","Cd86","Arg1",'Mrc1',"H2-Ab1","H2","Cd74")
   features =intersect(features,rownames(CARD_obj@refined_expression))
   p7 <- CARD.visualize.gene(
      spatial_expression = CARD_obj@refined_expression,
      spatial_location = location_imputation,
      gene.visualize = features,
      colors = NULL,
      NumCols = 6)
   print(p7)
   
   
   lapply(immune.marker,function(x){
      p=DotPlot(obj,features = x)
      data= p@
         scale_color_gradientn(colours = rcolors::rcolors$NCV_jet)+ggtitle(names(x))
   })
   
   features = c("Ly6g","Itgam","Itgax","Sigel","Cd14","Adgre1" ,"Cd163","Cd86","Arg1",'Mrc1',"H2-Ab1","H2","Cd74")
   features =intersect(features,rownames(CARD_obj@spatial_countMat))
   p8 <- CARD.visualize.gene(
      spatial_expression = CARD_obj@spatial_countMat,
      spatial_location = CARD_obj@spatial_location,
      gene.visualize = features,
      colors = NULL,
      NumCols = 6)
   print(p8)
   dev.off()
   return(CARD_obj)
   
})

 # #Extension of CARD for single cell resolution mapping 
 # scMapping = CARD_SCMapping(CARD_obj,shapeSpot="Square",numCell=20,ncore=10)
 # 
 # library(SingleCellExperiment)
 # MapCellCords = as.data.frame(colData(scMapping))
 # count_SC = assays(scMapping)$counts
 # df = MapCellCords
 # p10 = ggplot(df, aes(x = x, y = y, colour = CT)) + 
 #    geom_point(size = 3.0) +
 #    scale_colour_manual(values =  colors) +
 #    #facet_wrap(~Method,ncol = 2,nrow = 3) + 
 #    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
 #          panel.background = element_rect(colour = "white", fill="white"),
 #          plot.background = element_rect(colour = "white", fill="white"),
 #          legend.position="bottom",
 #          panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
 #          axis.text =element_blank(),
 #          axis.ticks =element_blank(),
 #          axis.title =element_blank(),
 #          legend.title=element_text(size = 13,face="bold"),
 #          legend.text=element_text(size = 12),
 #          legend.key = element_rect(colour = "transparent", fill = "white"),
 #          legend.key.size = unit(0.45, 'cm'),
 #          strip.text = element_text(size = 15,face="bold"))+
 #    guides(color=guide_legend(title="Cell Type"))
 # print(p10)
 # 

#cell marker variation
BiocManager::install(c('rhdf5','motifmatchr','chromVAR','Rsamtools','ComplexHeatmap'))
devtools::install_github("GreenleafLab/ArchR")
subtypecolor =ArchRPalettes$ironMan
names(subtypecolor)=0:14
p = DoHeatmap(immune.combined,features = unique(unlist(immune.markers)))
data = p$data
data = split(data,data$Feature)
cellmeta = immune.combined@meta.data
cellmeta$Cell = names(immune.combined@active.ident)

plt = lapply(1:length(data),function(i){
   df = data[[i]]
   df = left_join(cellmeta,df,"Cell")
   ggplot(df,aes(x=distance,y=Expression))+geom_point(data=df,aes(fill=Identity),shape=21)+
      geom_jitter(data=df,aes(fill=Identity),shape=21)+scale_fill_manual(values = subtypecolor)+
      geom_smooth(data=df,aes(color=Group))+ggtitle(names(data)[i])+theme_classic()+theme(legend.position = "bottom")
   
})


plt2 = lapply(1:length(data),function(i){
   df = data[[i]]
   df = left_join(cellmeta,df,"Cell")
   ggplot(df,aes(x=distance,y=Expression))+
      # geom_point(data=df,aes(fill=Identity),shape=21,size=0.1)+
      # geom_jitter(data=df,aes(fill=Identity),shape=21,size=0.1)+scale_fill_manual(values = subtypecolor)+
      geom_smooth(data=df,aes(color=Group))+ggtitle(names(data)[i])+theme_classic()+theme(legend.position = "top")
   
})
pdf("Marker_expression_distance.pdf",width=5,height = 4)
plt
dev.off()

pdf("Marker_expression_distance_nojitter.pdf",width=5,height = 3)
plt2
dev.off()



plt3 = lapply(1:length(data),function(i){
   df = data[[i]]
   df = left_join(cellmeta,df,"Cell")
   ggplot(df[df$Expression>0,],aes(x=distance,y=Expression))+geom_point(shape=21,aes(fill=Identity),)+
      geom_jitter(aes(fill=Identity),shape=21)+scale_fill_manual(values = subtypecolor)+
      geom_smooth(aes(color=Group))+ggtitle(names(data)[i])+theme_classic()+theme(legend.position = "bottom")
   
})

#====== Fig1 Recluster Global cells =====================

#read meta.data
setwd("~/LJX/sc_in")
obs = readr::read_csv("obs/obs_all_20230815.csv")
setwd("~/LJX/sc_in/")
n9 = Read10X("./N9/")
r5 = Read10X("./R5/")

obj.list = list('Nonresponder'=n9,'Responder'=r5)
# obj.list = lapply(obj.list,function(x){
#    obj = Read10X(x)
# })

n9.cd = obj.list$Nonresponder[,obs$barcodeID[obs$sample=="N9"]]
colnames(n9.cd)=paste0(colnames(n9.cd),"-1")

r5.cd = obj.list$Responder[,obs$barcodeID[obs$sample=="R5"]]
colnames(r5.cd)=paste0(colnames(r5.cd),"-2")


obs = data.frame("CellID"=obs$barcodeID,"sample"=obs$sample,"celltype"=obs$Cell_type,"cluster"=obs$CLUSTERS)
obs$CellID = ifelse(obs$sample=="N9", paste0(obs$CellID,"-1"),paste0(obs$CellID,"-2"))
obs$celltype =factor(obs$celltype,levels = unique(obs$celltype))

rownames(obs)=obs$CellID

cd = cbind(n9.cd,r5.cd)

all.equal(colnames(cd),obs$CellID)
cd = cd[,obs$CellID]

adta.obj = CreateSeuratObject(counts = cd,meta.data = obs)
#-----use nmf cluster-------
BiocManager::install('Biobase')
install.packages('NMF')

library(Biobase)
library(NMF)
# ## 参数测试
# data("esGolub")
# #esGolub <- esGolub[1:500,]
# t1 <- nmf(esGolub, 3, method = "brunet", seed = 219)
# runtime(t1)   # elapsed: 1.239 
# t2 <- nmf(esGolub, 3, method = "lee", seed = 219)
# runtime(t2)   # elapsed: 1.471 
# t3 <- nmf(esGolub, 3, method = "snmf/r", seed = 219)
# runtime(t3)   # elapsed: 0.995 
# t4 <- nmf(esGolub, 3, method = "brunet", seed = 'nndsvd')
# runtime(t4)   # elapsed: 3.156
# t5 <- nmf(esGolub, 3, method = "brunet", nrun = 100)
# runtime(t5)   # elapsed: 2.363
library(Seurat)
library(tidyverse)
library(NMF)
rm(list = ls())

## 创建seurat对象
# pbmc <- Read10X_h5("pbmc.h5")
# pbmc <- CreateSeuratObject(pbmc, project = "pbmc", min.cells = 3, min.features = 200)

#---use random 2000 cells test ----------

pbmc <- subset(adta.obj,cells = names(adta.obj@active.ident)[sample(2000)])
pbmc$percent.mt <- PercentageFeatureSet(pbmc, pattern = "^mt-")
pbmc <- subset(pbmc, percent.mt<20)
pbmc <- NormalizeData(pbmc) %>% FindVariableFeatures() %>% ScaleData(do.center = F)

## 使用pca的分解结果降维聚类
pbmc <- RunPCA(pbmc)
set.seed(219)
pbmc.pca <- RunUMAP(pbmc, dims = 1:20) %>% FindNeighbors(dims = 1:20) %>% FindClusters()

## 结果可视化
features= c("EGFP","mCherry-CAR","Cd3d","Cd4","Cd8a",
                 "Sell","Gzma",
                 "Gimap3","Il7r","Cd44","Gzmb","Prf1",
                 "Klrb1c","Nkg7","Mki67",
                 "Dcn","Col1a1","Col1a2","Col3a1",
                 "Lyz2","H2-Aa",'H2-Ab1',"EGFP",
                 "Cd19","Cd79a","Cd79b","Ighd",
                 "Foxp3","Il2ra",
                 'Epcam','Wfdc2','Slc12a2',
                 "Siglech")

features = c("EGFP","Cd8b1",
             "Ccr7","Fas",
             "Klrb1c","Nkg7","Mki67","Gzmb","Prf1",
             "Dcn","Col1a1","Col1a2","Col3a1",
             'Epcam','Wfdc2','Slc12a2',
             "Cd19","Cd79a","Cd79b","Ighd","Prg2",
             "Lyz2","H2-Aa",'H2-Ab1',"Cd14",
             "Adgre1","Mrc1","Itgax","Siglech"
             )

p <- DimPlot(pbmc.pca, label = T) + ggsci::scale_color_igv()
ggsave("pbmc_pca.png", p, width = 9, height = 6)
p <- FeaturePlot(pbmc.pca, features = features ,ncol = 8)
ggsave("pbmc_pca_markers.png", p, width = 12, height = 8)


## 高变基因表达矩阵的分解
# pbmc大体可分成T，B，NK，CD14+Mono，CD16+Mono，DC，Platelet等类型，考虑冗余后设置rank=10
vm <- pbmc@assays$RNA@scale.data
res <- nmf(vm, 10, method = "snmf/r", seed = 'nndsvd') 
runtime(res)
#    用户     系统     流逝 
#1063.147   78.019 1139.831 

## 分解结果返回suerat对象
pbmc@reductions$nmf <- pbmc@reductions$pca
pbmc@reductions$nmf@cell.embeddings <- t(coef(res))    
pbmc@reductions$nmf@feature.loadings <- basis(res)  

## 使用nmf的分解结果降维聚类
set.seed(219)
pbmc.nmf <- RunUMAP(pbmc, reduction = 'nmf', dims = 1:20) %>% 
   FindNeighbors(reduction = 'nmf', dims = 1:20) %>% FindClusters()

## 结果可视化  
p <- DimPlot(pbmc.nmf, label = T) + ggsci::scale_color_igv()
ggsave("pbmc_nmf.png", p, width = 9, height = 6)
p <- FeaturePlot(pbmc.nmf, features = features, ncol = 8)

ggsave("pbmc_nmf_markers.png", p, width = 12, height = 8)

# #------use harmony -------------
# install.packages('harmony')
# library(harmony)
# 
# immune.combined <- RunHarmony(immune.combined, "seurat_clusters")
# immune.combined <- RunUMAP(immune.combined, reduction = "harmony")


#old
adta.obj@active.ident = obs$celltype

markers = FindAllMarkers(adta.obj)

adta.obj = SCTransform(adta.obj)
adta.obj[["percent.mt"]] <- PercentageFeatureSet(adta.obj, pattern = "^mt-")
adta.obj<- subset(adta.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 5)

all.genes <- rownames(adta.obj)
adta.obj <- ScaleData(adta.obj, features = all.genes)
#    
# adta.allmarker = FindAllMarkers(adta.obj,only.pos = T,min.pct = 0.5)
# adta.allmarker = FindAllMarkers(adta.obj)
# adta.allmarkers = 

# features = adta.allmarker %>% group_by("celltype") %>%top_n

adta.obj <- FindVariableFeatures(adta.obj, x.low.cutoff = 0.0125, y.cutoff = 0.25, do.plot=FALSE)
adta.obj <- ScaleData(adta.obj, features =c("EGFP","mCherry-CAR",VariableFeatures(adta.obj)))

myfindcluster= function(obj,features,npcs,res,dim){
   pbmc <- obj
   pbmc <- RunPCA(pbmc,features = features,npcs = npcs)
   pbmc <- FindNeighbors(pbmc,dim=1:dim)
   pbmc <- FindClusters(pbmc,resolution = res)
   pbmc <- RunUMAP(pbmc,dim=1:dim)
   pbmc <- RunTSNE(pbmc)
}


adta.obj.immune = subset(sdta.obj,cells = obs$cellID[obs$celltype!="Tumor"])
adta.obj.immune = myfindcluster(adta.obj.immune,features =c("mCherry-CAR",features),npcs = 50,res = 1,dim = 15)


features =c(intersect(VariableFeatures(adta.obj),immunemarker$gene),"EGFP","mCherry-CAR") %>%unique()
features = VariableFeatures(adta.obj)

adta.obj <- myfindcluster(adta.obj,features = features,npcs = 50,res = 1,dim = 30)

# adta.obj.tumor = subset(adta.obj,cells = obs$CellID[obs$celltype=="Tumor"])

# adta.obj.tumor = FindVariableFeatures(adta.obj.tumor)

# adta.obj.tumor = myfindcluster(adta.obj.tumor,features = VariableFeatures(adta.obj.tumor),npcs = 50,res = 0.1,dim = 15)
# adta.obj<- merge (adta.obj.immune, y = adta.obj.tumor, add.cell.ids = c ("normal", "tumor"), project = "merged")

adta.obj <-FindVariableFeatures(adta.obj)


adta.obj.markers = FindAllMarkers(adta.obj)
#----figures--------
pdf("tc_features.pdf",width=10,height=8)
FeaturePlot(adta.obj,features = c("Cd3d","mCherry-CAR","Cd4","Cd8a"),label = T)
dev.off()

pdf('testmarkers.pdf',width = 10,height = 6)
DotPlot(adta.obj,features = c(
   "Cd3d","Cd4","Cd8a",
   "Sell","Gzma",
   "Gimap3","Il7r","Cd44","Gzmb","Prf1",
   "Klrb1c","Nkg7","Mki67",
   "Dcn","Col1a1","Col1a2","Col3a1",
   "Lyz2","H2-Aa",'H2-Ab1',"EGFP",
   "Cd19","Cd79a","Cd79b","Ighd",
   "Foxp3","Il2ra",
   'Epcam','Wfdc2','Slc12a2',
   "Siglech"
    )
        # group.by = "celltype"
        )+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

# DotPlot(adta.obj,features = c(
#                               "Cd3d","Cd4","Cd8a",
#                               "Sell","Gzma",
#                               "Gimap3","Il7r","Cd44","Gzmb","Prf1",
#                               "Klrb1c","Nkg7","Mki67",
#                               "Dcn","Col1a1","Col1a2","Col3a1",
#                               "EGFP","Lyz2",
#                               "Cd19","Cd79a","Cd79b","Ighd",
#                               "Foxp3","Il2ra",
#                               'Epcam','Lgals3','Wfdc2','Slc12a2',
#                               "Siglech"
#                               
#                               
#                               ),
#         group.by = "celltype"
# )+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()

#check previous annotation
pdf('test_umap_split.pdf',width=10,height = 5)
p1=DimPlot(adta.obj,group.by = "celltype",split.by = "sample",label = T)

p2=DimPlot(adta.obj,label = T,split.by = "sample")

p1 
p2


dev.off()

markers =FindAllMarkers(adta.obj,only.pos = T,min.pct =0.5)

markers[markers$cl]

active.ident = adta.obj@active.ident
celltypename =ifelse(active.ident%in%c(0,1,))
celltypename = ifelse(active.ident=="14","Tumor",as.character(obs$celltype))


# markers_celltype =FindAllMarkers(adta.obj,group.by="celltype",only.pos = T,min.pct = 0.5)
#new annotation
rename = c(
   '0'='DPT',
   '1'='DPT',
   '2'='DPT',
   '3'='CD4T',
   '4'='CD4T',
   '5'='CD4T',
   '6'='CD4T',
   '7'='DPT',
   '8'='CD8T',
   '9'="DNT",
   '10'="NKT",
   '11'="Cycling T",
   '13'='CAFs',
   '12'="Myeloid",
   "14"="Tumor",
   "15"="B Cell",
   "16"="DNT",
   "17"="Treg",
   '18'="Cycling T",
   '19'='Epithelium',
   '20'='CAFs',
   '21'='pDC'
   )

adta.obj = RenameIdents(adta.obj,rename)

pdf('test_umap_split.pdf',width=12,height = 5)

p1=DimPlot(adta.obj,split.by = "sample",label = T,repel = T)
df =p1$data 
length(levels(adta.obj@active.ident))
cellcolors = ArchR::ArchRPalettes$circus


names(cellcolors)
colors = cellcolors[1:length(levels(adta.obj@active.ident))]
names(colors)=levels(adta.obj@active.ident)
# df$ident = factor(df$ident,levels=names(colors))

mylabelcellnum=function(p){
   totalcd4 = nrow(p$data)
   cellnum =table(p$data$ident)%>%as.data.frame()
   legendlabel = paste0(cellnum$Var1,"(",cellnum$Freq,")")
   df=p$data
   df$ident =factor(df$ident,
                    levels = cellnum$Var1,
                    labels = legendlabel)
   return(df)
}
df=mylabelcellnum(p1)
df$sample = factor(df$sample,levels=c("N9","R5"),labels = c("Non-Response","Response"))

pdf("test_umap.pdf",width=6,height = 4)
names(colors)=levels(df$ident)
p1 = ggplot(df,aes(x=UMAP_1,y=UMAP_2))+
   geom_point(aes(color=ident),size=0.3)+scale_color_manual(values = colors)+
   # facet_grid(.~sample)+
   cowplot::theme_map()

LabelClusters(p1,id = "ident",box = T,color="white")+theme(text = element_text(size = 14))+NoLegend()

dev.off()

cellratio =table(df$ident,df$sample)%>%as.data.frame()
samplenum=table(df$sample)
cellratio$Var2 =factor(cellratio$Var2,levels=c("N9","R5"),labels = c(paste0("Non-response\n(",samplenum[[1]],")"),
                                                                          paste0("Response\n(",samplenum[[2]],")") ) )
pdf("./test_cellratio.pdf",width =6,height = 4)
ggplot(cellratio,aes(x=Freq,y=Var2,fill=Var1))+
   geom_bar(stat = "identity",position = "fill")+
   scale_fill_manual(values = colors,name="")+theme_pubr(border = T,legend = "right")+xlab("% Total")+ylab("")+
   coord_flip()
dev.off()

pdf('testmarkers.pdf',width = 10,height = 6)
DotPlot(adta.obj,features = c(
   "Cd3d","Cd4","Cd8a",
   "Sell","Gzma",
   "Gimap3","Il7r","Cd44","Gzmb","Prf1",
   "Klrb1c","Nkg7","Mki67",
   "Dcn","Col1a1","Col1a2","Col3a1",
   "Lyz2","H2-Aa",'H2-Ab1',"EGFP",
   "Cd19","Cd79a","Cd79b","Ighd",
   "Foxp3","Il2ra",
   'Epcam','Wfdc2','Slc12a2',
   "Siglech"
)
# group.by = "celltype"
)+theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

dev.off()
   # ---=ComplexHeamap for celltype top markers ----
   
   markers = FindAllMarkers(adta.obj,only.pos = T,min.pct = 0.2)
   markers$pct.fc = markers$pct.1/markers$pct.2
   topgene = markers[markers$pct.1>5*markers$pct.2,]
   table(topgene$cluster)
   topgene = markers %>% group_by(cluster) %>% top_n(pct.fc,n=10)
   p = DoHeatmap(adta.obj,features = topgene$gene,group.colors = colors,)
   mat = Matrix(data=p$data$Expression,
                nrow = length(levels(p$data$Feature)),
                ncol = length(levels(p$data$Cell)),dimnames = list(levels(p$data$Feature),levels(p$data$Cell)))
   
   mat = as.matrix(mat)
   all.equal(colnames(mat),names(adta.obj@active.ident))
   #cellratio
   cellratio = table(adta.obj@active.ident)
   
   mat = mat[,names(adta.obj@active.ident)]
   library(ComplexHeatmap)
   annotation = data.frame(adta.obj@active.ident)
   colnames(annotation)= c("Cell Type")
   pdf("test_complexheatmap_markers_celltype.pdf",width = 10,height = 10)
   is.na(mat)<-0
   p=ComplexHeatmap::pheatmap(mat,show_rownames = F,
                              show_colnames = F,scale="row",
                              color = rev(rcolors$RdBu),
                              treeheight_row = 0,cluster_cols = F,
                              annotation_col = annotation,
                              annotation_colors = list("Cell Type"=colors),
                              column_split = adta.obj@active.ident)
   A=p@matrix
   
   genes = c(
      "Cd3d","Cd4","Cd8a",
      "Sell","Gzma",
      "Gimap3","Il7r","Cd44","Gzmb","Prf1",
      "Klrb1c","Nkg7","Mki67",
      "Dcn","Col1a1","Col1a2","Col3a1",
      "Lyz2","H2-Aa",'H2-Ab1',"EGFP",
      "Cd19","Cd79a","Cd79b","Ighd",
      "Foxp3","Il2ra",
      'Epcam','Wfdc2','Slc12a2',
      "Siglech")
   genes = data.frame(genes=genes)
   
   p=p+rowAnnotation(link = anno_mark(at = which(rownames(A) %in% genes$genes), 
                                    labels = genes$genes, 
                                    labels_gp = gpar(fontsize = 10)))
   
   print(p)
   dev.off()
   
   png("Heatmap_features.png",width = 2000,height = 1600,res = 300)
   DoHeatmap(adta.obj,features = topgene$gene,group.colors = colors)
   p+scale_fill_gradientn(colors= rev(rcolors::rcolors$RdBu))
   
   dev.off()

   
   
#-----RunCCA & RunsPCA for ST data--------
   
   he_fig.n9 ="~/LJX/st/n9/05.AllheStat/allhe/he_roi_small.png"
   he_fig.r5 ="~/LJX/st/r5/05.AllheStat/allhe/he_roi_small.png"
   
   levelpath = dir(n9_path)
   #----- create reference scRNA data -----------
   reference$celltype = reference@active.ident
   immune.combined =reference
   
   reference = subset(immune.combined,cells =names(immune.combined@active.ident)[immune.combined@active.ident%in%c('EpC','Endo','Fibro','Tumor')])

   DefaultAssay(reference)<-'RNA'
   reference = SCTransform(reference)
   
   
   #----------create query ST seurat object -------------------
   st.n9=mySeurat(n9_path,levelpath = levelpath[3],
                  hefigpath = he_fig.n9,projname = "NR",
                  # features = c('EGFP','mCherry-CAR','Cd19'),
                  min.features = 50,nfeatures = 2000)
   
   
   st.r5=mySeurat(r5_path,levelpath = levelpath[3],
                  hefigpath = he_fig.r5,projname = "R",
                  # features = c('EGFP','mCherry-CAR','Cd19'),
                  min.features = 50,nfeatures = 2000)
   
   st.list =list("NR"=st.n9,
                 "R"=st.r5)
   
   cca.n9 = RunCCA(object1 = , object2 = pbmc2)
   
   myAddspatial = function(pbmc){
     pbmc@reductions$spatial <- pbmc@reductions$umap
     spatial_embedding = data.frame(pbmc@meta.data[,c("x","y")])
     colnames(spatial_embedding)=c("SPATIAL_1","SPATIAL_2")
     pbmc@reductions$spatial@cell.embeddings <- as.matrix(spatial_embedding)
     pbmc@reductions$spatial@key="SPATIAL_"
     return(pbmc)
   }
   
   
   st.list = lapply(1:length(st.list),function(x)myAddspatial(st.list[[x]]))
   # pbmc@reductions$spatial@feature.loadings <- basis(res)
   
  #-----------CCA mapping -----------------
   myCCAmap=function(obj.query,obj.reference,plot=FALSE){
     pancreas.query <- obj.query
     pancreas.anchors <- FindTransferAnchors(reference =obj.reference, 
                                             query = obj.query,
                                             dims = 1:30, reference.reduction = "pca")
     predictions <- TransferData(anchorset = pancreas.anchors, refdata = obj.reference$celltype,
                                 dims = 1:30)
     pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
     
     if(plot==TRUE){
       p1 <- DimPlot(obj.inference, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
                     repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
       p2 <- DimPlot(pancreas.query, reduction = "umap", group.by = "predicted.celltype", label = TRUE,
                     label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
       p1 + p2
       print(p1+p2)
       
     }else{
       return(pancreas.query)}
   }

   # reference$celltype = reference@active.ident
   
   reference$celltype =factor(reference$celltype,
                              levels=unique(reference$celltype))
   reference.list <- SplitObject (reference, split.by = "sample")
   st.n9 = myCCAmap(obj.query = st.n9,obj.reference = reference.list$Nonresponder)
   st.r5 = myCCAmap(obj.query = st.r5,obj.reference = reference.list$Responder)
   plotCCAmap =function(obj.inference,obj.query,group_by='predicted.id',reduction='umap',cols=NULL){
     pancreas.query=obj.query
     if(is.null(cols)){
       p1 <- DimPlot(obj.inference, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
                     repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
       p2 <- DimPlot(pancreas.query, reduction = reduction, group.by = group_by, label = TRUE,
                     label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
     }else{
       p1 <- DimPlot(obj.inference, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
                     repel = TRUE,cols=cols) + NoLegend() + ggtitle("Reference annotations")
       p2 <- DimPlot(pancreas.query, reduction = reduction, group.by = group_by, label = TRUE,
                     label.size = 3, repel = TRUE,cols = cols) + NoLegend() + ggtitle("Query transferred labels")
       
     }
          return(p1 + p2)
   }
   
   plotCCAmap(obj.inference = reference.list$Nonresponder,obj.query = st.n9)
   plotCCAmap(obj.inference = reference.list$Responder,obj.query = st.r5,reduction = 'spatial',cols = ctcolors)
   
   st.list = lapply(1:length(st.list),function(i){
     obj.query = st.list[[i]]
     obj.query =myCCAmap(obj.query = obj.query,
                         obj.reference = immune.combined)
     return(obj.query)
   })
   
   names(st.list)=c('Nonresponder','Responder')
  # st.n9 = myCCAmap(obj.query = st.n9,
  #                  obj.reference = immune.combined)
  # 
  # st.r5 = myCCAmap(obj.query = st.r5,
  #                  obj.reference = immune.combined) 
   myFeaturePlot = edit(myFeaturePlot)
   themes = cowplot::theme_map()
   plt =lapply(st.list,function(x){
     obj =x
     p = FeaturePlot(x,features = c("EGFP","mCherry-CAR",
                                    grep("score",colnames(x@meta.data),value = T)),
                     reduction = "spatial",
                     # xlim = c(0,1000),ylim = c(0,1000),
                     # geom_theme = themes,
                     cols = rcolors::rcolors$cmocean_ice,pt.size = 0.3,)&xlim(1,1000)&ylim(1,1000)
       # xlim(0,1000)+ylim(0,1000)

     return(p)   })
   
   ggsave("spatialplot_scrna_cca.png",CombinePlots(plt),width = 35,height = 20)
   
   
#--------use classical markers------------------
#ST.NR
   DimPlot(st.n9,reduction = 'spatial',
           cols = rcolors::rcolors$GHRSST_anomaly)
   st.nr.name = c('0'='TAM',#Lyz2,Ctss,H2-D1,Cd74
                  '1'='T',#Cox8a,B2m
                  '2'='Erythoid',#Hbb-bs,Hba-a1,Hba-a2
                  '3'='',
                  '4'='Dead Epi',#Filip1l induce cell apoti
                  '6'='TAM/Fib',#C1qa,Mrc1,Apoe
                  '17'='Neutrophi',#S100a8,S100a9
                  '23'='Inflammation', #Cxcl2,Cxcl3,Ptgs2
                  '25'='Erythoid',
                  '33'='B',#Ighm,Igkc
                  '35'='CAF'
                  )
   
   Celltype = list('SPP1+ TAM'=c('0'),
                   'Tumor'=c('1'),#EGFP,B2m
                   'Vessel'=c('2'),#Pecam1+Vwf+Hba-a1+
                   'Pro-Vessel'=c('25'),#Igfbp7,Col4a1
                   'Endo'=c('10'),#Pecam1+, Igfbp7,Col4a1,mesenchymal
                   'Dead Tumor'=c('4'),#?Filip1l,epithelial-mesenchymal,Cell Apoptosis
                   'Pro-Fibo'=c('18'),#Igfbp7,Col4a1
                   'MDSC/Fib'=c('5'),#Sparc,Col1a1,Col1a2,Col3a1,Ptx3
                   'LAM/Fib'=c('6'),#Col1a2,Col1a1,C1qa,Mrc1,Apoe,Saa3
                   'Neutro'=c('17'),#S100a9,S100a8
                   'Erythoid'=c('9','15'),
                   'Inflammation'=c('23'),#Cxcl2,Cxcl3
                   'B'=c('33'),
                   'CAF'=c('35'),#Serpinb2,Ptx3,Col1a1C lipid-asscoiated
                   'microbe mt-related'=c(
                     '3',#mt-Atp6,'mt-Co1' disease,Cytochrome C Oxidase I
                     
                     '8','13','14','16','30','34',#mt-Rnr2 
                                          '11','19','20','21','22','27','29','31','32'),#Gm10800
                   'Ribosom-related'=c('12',
                                       '24' #RBM22
                                       ),
                   'Unknown'=c('26',#Dnah7b+
                               '28',#Gm19951
                               '36' ,#Gm10719
                               '7',#CT010467.1
                               '13'#Camk1d
                               )
                   )
   
   
  niche.n9 = reshape2::melt(Celltype) 
  niche.n9$L1 = factor(niche.n9$L1,levels=unique(niche.n9$L1))
  niche.n9$value=factor(niche.n9$value,levels = levels(st.n9@active.ident))
  niche = niche.n9$L1
  names(niche)=niche.n9$value
  st.n9@active.ident = st.n9$seurat_clusters
  st.n9 = RenameIdents(st.n9,niche)
   
  unknown = names(st.n9@active.ident)[st.n9@active.ident%in%c('Unknown',
                                                            'Ribosom-related',
                                                            'microbe mt-related')]
  known = names(st.n9@active.ident)[!st.n9@active.ident%in%c('Unknown',
                                                            'Ribosom-related',
                                                            'microbe mt-related')]
  known =subset(st.n9,cells = known)

  unknown = subset(st.n9,cells =unknown)
  
  unknown[["percent.mt"]] <- PercentageFeatureSet(unknown, pattern = "^mt-")
  
  unknown <- subset(unknown, 
                    subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & 
                      percent.mt < 5)
  
  # unknown = subset(st.n9,cells= )
  
   DotPlot(st.n9,features = c(
     'EGFP',
     'Cd3e','Cd3d','Cd3g','Cd5','Cd7',
     'Hbb-bt','Hba-a1','Hba-a2','Hbb-bs',
                              'Cxcl2','Cxcl3',
                              'C1qa','C1qc','Mrc1','Apoe','Spp1','Tgfb1',
                              'Ighm','Igkc',
                              'S100a8','S100a9',
                              'Col1a1','Col1a2','Dcn',
                              'Igfbp7','Col4a1','Gm10800'),cluster.idents = T)+coord_flip()
   
#=======Fig2 Global T cell reclassification and comparison with CAR-T========================
library(AUCell)
library(clusterProfiler)
setwd("~/LJX/sc_in/")
dir.create("TC_subcluster")
setwd("TC_subcluster")

#---previous read ref_celltype from CellMarker 2.0 
ref_celltype=readxl::read_xlsx("~/LJX/mm/Cell_marker_Mouse.xlsx")

head(ref_celltype)

ref_celltype =data.frame("term"=ref_celltype$cell_name,"gene"=ref_celltype$Symbol)

tc_celltype = ref_celltype[grep("T ",ref_celltype$term),]

tc_meta = sc_meta[grep("T",sc_meta$cellType),]
tc_meta = tc_meta[grep("Tumor",tc_meta$cellType,invert = T ),]
tc_count = sc_count[,tc_meta$cellID]

set.seed(12345)

#global T
tccells = names(adta.obj@active.ident)[adta.obj@active.ident%in%c("CD4T","DPT", "DNT",  "Cycling T", "Treg","CD8T","NKT")]
adta.obj@active.ident[tccells] %>% unique()
tc_obj = subset(adta.obj,cells= tccells)
saveRDS(adta.obj,"~/LJX/sc_in/test_seurat_all.RDS")

rm(adta.obj)

gc()
expr_matrix =tc_obj@assays$RNA@counts
gene_annotation=data.frame("gene_id"=rownames(expr_matrix),
                           "gene_short_name"=rownames(expr_matrix))

rownames(gene_annotation)=gene_annotation$gene_id
sample_sheet = tc_obj@meta.data
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

# pd <- new("AnnotatedDataFrame", data = pd)
# fd <- new("AnnotatedDataFrame", data = fd)

cds <- newCellDataSet(expr_matrix, #SparseMatrix format
                      phenoData = pd, featureData = fd)

# rpc_matrix <- relative2abs(cds, method = "num_genes")

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
gc()

cds <- detectGenes(cds, min_expr = 0.1)
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) 
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
diff_test_res <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters",
                                      cores = 10)

ordering_genes <- row.names (subset(diff_test_res, qval < 0.05))
cds <- setOrderingFilter(cds, ordering_genes)
cds = reduceDimension(cds, max_components = 2,
                      method = 'DDRTree')

cds = order(cds)
#old
tc_obj = mypreprocess(count = tc_count,meta.data = tc_meta,res = 1,dim=30)

tc_markers= FindAllMarkers(tc_obj,only.pos = T,logfc.threshold = 0.1,min.pct = 0.3)
tc_markers_sub =tc_markers[tc_markers$pct.1>2*tc_markers$pct.2,]
# write.csv(celltype,"tc_enrichr_sub_cellmarker.csv")

celltype = myenrichr(markers = tc_markers_sub,ref_celltype = tc_celltype)
write.csv(celltype,"tc_enrichr_sub_cellmarker.csv")
tc_markers_sub =left_join(tc_markers_sub,ref_celltype,"gene")

# tc_subtypes = myenrichr(markers = tc_markers)

# tc_obj =subset(tc_obj,cells = names(tc_obj@active.ident)[grep("16",tc_obj@active.ident,invert = T)])
# tc_obj =subset(tc_obj,cells = names(tc_obj@active.ident)[grep("21",tc_obj@active.ident,invert = T)])
# tc_obj =subset(tc_obj,cells = names(tc_obj@active.ident)[grep("20",tc_obj@active.ident,invert = T)])
# 
# megatc = list(DNT=c("1","7","8","9","11","12","15","18"),
#    DPT=c("0","2","4","14","17"),
#    CD4T=c("3","6","13"),
#    CD8T=c("5","10","19"))
# megatc =reshape2::melt(megatc)
# clusters = as.data.frame(tc_obj@active.ident)
# head(megatc)
# colnames(megatc)=c("cluster","MegaTC")
# colnames(clusters)=c("cluster")
# clusters = left_join(clusters,megatc)
# clusters$CellID = names(tc_obj@active.ident)
# tc_obj =AddMetaData(tc_obj,clusters)
# saveRDS(tc_obj,"tc_mega.RDS")
# clusters = tc_obj@meta.data

tc_obj$cellType =factor(tc_obj$cellType,levels=c("DNT","DPT","CD4T Cell","CD8T Cell","T Cell MKI67+","NKT"))
clusters = data.frame("CellID"=names(tc_obj@active.ident),"cluster"=tc_obj@active.ident,"MegaTC" = tc_obj$cellType)
tc.list = split(clusters$CellID,clusters$MegaTC)
tc_sub.list = lapply(tc.list,function(x){
   obj.sub = subset(tc_obj,cells=x)
   count =  obj.sub@assays$RNA@counts
   meta = obj.sub@meta.data
   obj.sub = mypreprocess(count =count,meta.data = meta,res = 1)
})

tc_obj_coord =cbind(tc_obj@reductions$umap@cell.embeddings,tc_obj@reductions$tsne@cell.embeddings)
rm(tc_obj)
gc()

tc_submarkers =lapply(tc_sub.list,function(x){
   obj = x
   markers = FindAllMarkers(obj,only.pos = T,logfc.threshold = 0.1)
   markers = markers[markers$pct.1>markers$pct.2,]
   return(markers)
})

tc_submarkers =data.table::rbindlist(tc_submarkers,
                                        use.names = T,
                                        fill = T,idcol = "MegaTC")
write.csv(tc_submarkers,"tc_subcluster_markers.fine.csv")

tc_umap.list = lapply(1:length(tc_sub.list),function(x){
   DimPlot(tc_sub.list[[x]],label=T)+ggtitle(names(tc_sub.list)[x])
})

pdf("tc_sub_umap.pdf",width = 15,height = 8)

CombinePlots(plots = tc_umap.list,ncol = 3,legend ='none' )
dev.off()

geneset = mm.genes[grep("T ",mm.genes$cell_name),]
geneset = geneset[grep("ST-HSC",geneset$cell_name,invert = T),]
geneset =geneset[,c("Symbol","cell_name")]
geneset =na.omit(geneset)
geneset =rbind(geneset,data.frame("Symbol"=c('Cd3e','Il22','Ahr'),"cell_name"=rep("Th22",3)))
geneset =split(geneset$Symbol,geneset$cell_name)
# geneset0 =list(geneset,"Th22"=c('Cd3e','Il22','Ahr'))
ref_celltype = reshape2::melt(geneset)
ref_celltype = data.frame("term"=ref_celltype$L1,"gene"=ref_celltype$value)

enrichr.celltypelist =lapply(1:length(tc_sub.list),function(x){
   res = myenrichrCelltype(obj=tc_sub.list[[x]],logfc = 0.1,
                           ref_celltype = ref_celltype,
                           tc.mega = names(tc_sub.list)[x])
   write.csv(res,paste0("CellMarker_enrichr_",names(tc_sub.list)[x],".csv"))
   return(res)
})

names(enrichr.celltypelist)=names(tc_sub.list)
enrichr.celltype =data.table::rbindlist(enrichr.celltypelist,use.names = T,fill = T,idcol = "MegaTC")


auc.list=lapply(tc_sub.list,function(x){
   obj = x
   expr=as.matrix(obj@assays$RNA@counts)
   tc_auc = myAUCell(exprMatrix =expr ,geneSets =geneset ,plot = T)
   dev.off()
   dev.new()
})


#----CAR-T variation---------------
tc_sub.list =lapply(tc_sub.list,function(x){
   mcherry=as.numeric(x@assays$RNA@counts["mCherry-CAR",])
   car=ifelse(mcherry >0,"CAR-T","Normal-T")
   x = AddMetaData(x,car,"CAR-T")
   return(x)
})

cellratio.list = lapply(tc_sub.list,function(x){
   table(x$sampleInfo,x$`CAR.T`) %>% as.matrix() %>% as.data.frame()
})
cellratio =data.table::rbindlist(cellratio.list,use.names = T,fill = T,idcol = "MegaTC")
cellratio$MegaTC = factor(cellratio$MegaTC,levels = names(cellratio.list))
devtools::install_github("XiaoLuo-boy/ggheatmap")
library(ggheatmap)
library(aplot)
library(ggrepel)
cellratio[cellratio$Var2=="CAR-T",]$Freq =100*cellratio[cellratio$Var2=="CAR-T",]$Freq/sum(cellratio[cellratio$Var2=="CAR-T",]$Freq)
cellratio[cellratio$Var2=="Normal-T",]$Freq =100*cellratio[cellratio$Var2=="Normal-T",]$Freq/sum(cellratio[cellratio$Var2=="Normal-T",]$Freq)

p=ggplot(cellratio,aes(x=MegaTC,y=Var1,fill=Freq),color="black")+
   geom_tile(color = "white")+scale_fill_viridis_c()+
   theme_minimal()+ # minimal theme
   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                    size = 10, hjust = 1),
         axis.text.y = element_text(
            size = 10)
   )+facet_wrap(.~Var2)+ coord_fixed()+ geom_text(aes(label =Freq),color="white")

#ggtree
ggheatmap(cellratio,aes(x=MegaTC,y=Var1,fill=Freq),color="black")+
   geom_tile(color = "white")+scale_fill_viridis_c()

+
   scale_fill_gradient2(low = "blue", high = "red", mid = "grey100", 
                        midpoint = 0.5, limit = c(0, 0.9), space = "Lab", 
                        name="Cell Count") +
   theme_minimal()+ # minimal theme
   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                    size = 10, hjust = 1),
         axis.text.y = element_text(
            size = 10)
   )+facet_wrap(.~Var2)+ coord_fixed() 
   # +geom_rect(xmin = 0.5, xmax = 3.5, ymin = 0.5, ymax = 3.5, col = "black", alpha = 0) 


pdf("tc_sub_umap.pdf",width = 8,height = 3)
lapply(1:length(tc_sub.list),function(x){
   p1= DimPlot(tc_sub.list[[x]],label = T,order = T,repel = T)+ggtitle(names(tc_sub.list)[x])
   p2 =DimPlot(tc_sub.list[[x]],label = T,group.by = "CAR.T",order = T,repel = T)
   p1+p2
})
dev.off()

#---- general monocle evolution -----------

tcs=names(immune.combined@active.ident)[immune.combined@active.ident%in%c("Tex", 
                                                                          "Tpex", 
                                                                          "Treg", 
                                                                          "NK/NKT",
                                                                          "CLT", 
                                                                          "Teff",  
                                                                          "Tem",  
                                                                          "Tcm",
                                                                          "Tn_CD69+","Tn")]
tc.obj = subset(immune.combined,cells=tcs)


library(parallel)
library(monocle)
n_cors = detectCores()

cl <- makeCluster(n_cors[1]-1)

setwd("~/LJX/sc_in/TC_subcluster/CD4T/")
sample_sheet=df #df is the data from processed Dimplot above
sample_sheet$ident = CD4T@active.ident
sample_sheet$group =factor(CD4T$sampleInfo,levels=c("N9","R5"),
                           labels = c("Non-response","Response"))
sample_sheet$CART_Abundance = sample_sheet$CART
sample_sheet$CART = factor(ifelse(sample_sheet$CART>0,"CAR-T","Normal-T"),
                           levels=c("CAR-T","Normal-T"))
head(sample_sheet)
str(sample_sheet)
expr_matrix=CD4T@assays$RNA@counts

gene_annotation=data.frame("gene_id"=rownames(expr_matrix),
                           "gene_short_name"=rownames(expr_matrix))

rownames(gene_annotation)=gene_annotation$gene_id
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)

cds <- newCellDataSet(expr_matrix, #SparseMatrix format
                      phenoData = pd, featureData = fd)

# rpc_matrix <- relative2abs(cds, method = "num_genes")

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
gc()
cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) 
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)

expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
# dispersionTable(cds)
diff_test_res <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = "~ident",
                                      cores = 10)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.05))
ordering_genes <- intersect(VariableFeatures(CD4T),ordering_genes)

cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
cds = reduceDimension(cds, max_components = 2,
                      method = 'DDRTree')

p1=plot_cell_trajectory(cds, color_by = "ident")+scale_color_brewer(palette = "Paired")
df = p1$data
p=ggplot(df,aes(x=data_dim_1,y=data_dim_2))+geom_point(aes(color=ident))+cowplot::theme_map()+scale_color_brewer(palette = "Paired")
pdf("trajectory_celltype.pdf",width=5,height = 4)
LabelClusters(p,id="ident",repel = T)
dev.off()

cds <- orderCells(cds)
#---- CellChat communication -------



#=========Fig3 Global B cell reclassification======================
setwd("~/LJX/sc_in/")
dir.create("BC_subcluster")
setwd("BC_subcluster")

geneset = mm.genes[grep("B ",mm.genes$cell_name),]
geneset =geneset[,c("Symbol","cell_name")]
geneset =na.omit(geneset)

bc_meta = sc_meta[grep("B",sc_meta$cellType),]
bc_count = sc_count[,bc_meta$cellID]

bc_obj = mypreprocess(count = bc_count,meta.data = bc_meta,res=1,dims)





#=========Fig4 Global myeloid cell reclassification================
setwd("~/LJX/sc_in/")
dir.create("Mye_subcluster")
setwd("Mye_subcluster")

mye_meta = sc_meta[grep("Mye",sc_meta$cellType),]
mye_count = sc_count[,mye_meta$cellID]

mye_obj = mypreprocess(count = mye_count,meta.data = mye_meta)

#========Fig5 Global fibroblast cell classification================
setwd("~/LJX/sc_in/")
dir.create("Mye_subcluster")
setwd("Mye_subcluster")

fib_meta = sc_meta[c(grep("CAF",sc_meta$cellType),
                     grep("fibroblast",sc_meta$cellType)),]

fib_count = sc_count[,fib_meta$cellID]

fib_obj = mypreprocess(count = fib_count,meta.data = mye_meta)



#========Fig6 Cell-cell interactions =============================




