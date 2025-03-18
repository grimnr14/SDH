#library(usdata)
library(tidycensus)
library(tidyverse)
library(data.table)
#rm(list = ls())
library(stringr)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(here)
library(usmap)
library(ggplot2)
library(psych)
library(dplyr)
#library(conflicted)
library(comorbidity)
library(table1)
library(liver)
library(fastICA)
library(rpart)
library(zipcodeR)
library(ipumsr)
#library(mice)
#library(sf)
#library(maptools)
#library(tigris)
#devtools::install_github('walkerke/tigris')
#options(tigris_use_cache=T)

makeIndex<-function(x,#x is a data.frame without NA values for features used in the factor analysis
                    geoid="geoid",#any coded value for geography level, which must be specified by variable name
                    state=NULL,#subsets values according to a string for lower case state name
                    cols=NULL,#optional vector of variable names to include in factor analysis, else the full data frame is used
                    nfactors=4,#number of factors solved for in setting the variable groupings
                    exclude=NULL,#a vector of variable names to exclude from final pca scoring of index
                    set.factor=NULL#a list of variable names to replace $factor_members, allowing a manual group assigment for scoring
){
  #Step 1: locate relevant features and load into a data.frame x that includes geography NAME and geoid----
  if(!is.null(state)){#if this is a state specific analysis, you may declare a value for state (full state name)
    x<-x[stringr::str_detect(tolower(x[,"NAME"]),state)==T,]
  }
  if(!is.null(cols)){
    x<-x[,c("NAME",geoid,cols)]
  }
  #Step 2: standardize inputs----
  x<-na.omit(x)#if missing ANYTHING, remove observation
  labs<-as.data.frame(cbind(x[,c("NAME",geoid)],rowid=rownames(x)))
  labs<-labs[order(labs$rowid),]
  x <- x %>% dplyr::mutate_at(vars(-c("NAME",geoid)), as.numeric)#ensure everything is numeric
  
  x<-x[,!names(x) %in% c(geoid,"NAME")]
  #x<- liver::transform(x, method = "zscore") #x is now normalized
  x<-liver::zscore(x)
  
  #Step 3: Explore and prune irrelevant, redundant features----
  corr_matrix <- cor(x)
  #p1<-corrplot::corrplot(corr_matrix, type = "upper")#get a bivariate corrplot for features
  
  res.pca <- FactoMineR::PCA(x,ncp=ncol(x),graph=FALSE)#conduct a pca
  # Eigenvalue
  eig.val <- data.frame(factoextra::get_eigenvalue(res.pca))
  t1<-round(eig.val,2)#store eigenvalues
  
  #Eigenvector
  eigenvectors<- res.pca$ind$coord
  t2<-eigenvectors#store eigenvectors for solution
  
  # Scree Plot, optional, using eigenvalue-1 criteria
  p2<-ggplot(data=eig.val,aes(y=eigenvalue,x=1:nrow(eig.val)))+
    geom_line(stat="identity")+
    geom_point(stat="identity")+
    geom_hline(yintercept=1,lty=2,col="darkred")+
    theme_minimal()+
    xlab("Dimensions")+
    ylab("Eigenvalue")
  
  # Plot the eigenvalues/variances against the number of dimensions
  var <- factoextra::get_pca_var(res.pca)
  
  # PCA plotting (Correlation circle)
  # - Positively correlated variables are grouped together.
  # - Negatively correlated variables are positioned on opposite sides
  # - Variables that are away from the origin are well represented on the factor map
  # - Color represents the contribution of variables
  # var$coord # correlation
  
  # var$cos2 # quality of representation
  #p3<-corrplot::corrplot(var$cos2, is.corr=FALSE)#stored correlation of eigenvectors
  
  # compute varimax-rotated principle components, the number of components can be specified.
  
  res.pca.rotated <- psych::principal(x, rotate="varimax",nfactors=nfactors,scores=T)
  #print(res.pca.rotated)#this is part of the output
  
  # Perform exploratory factor analysis
  # check each variable's loading against each factor
  fa_result <- psych::fa(x,nfactors = nfactors, rotate="varimax")
  #print(fa_result)#this is also part of the output
  
  #loadings<-fa_result$loadings#use loadings to create a fa diagram, except this does not store well and needs to be outputted
  
  #Check single solution
  fa_result2 <- psych::fa(x,nfactors = 1, rotate="varimax")
  #print(fa_result2)
  
  # Finalize the list of variables, remove the irrelevant indicators
  
  # Map Subcategories from fa loadings:
  ex<-fa_result$loadings
  outs<-NULL
  for(i in 1:nrow(ex)){
    outs<-c(outs,max(abs(ex[i,])))
  }
  maxs<-outs
  outs<-list()
  for(i in 1:nfactors){
    ex<-rownames(fa_result$loadings)[round(abs(fa_result$loadings[,i]),3)>0.3&
                                       !rownames(fa_result$loadings) %in% unlist(outs)&
                                       abs(fa_result$loadings[,i])==maxs]#we need some heuristic to decide loadings
    outs<-append(outs,list(ex))
  }
  factor_membership<-outs#store this for output
  # visualize the factor model
  #fa.diagram(fa_result$loadings)
  
  #Step 4: Generate the composite score using PCA----
  score <- x[,!names(x) %in% exclude]
  pca_result <- psych::principal(score)
  final_x <- as.data.frame(unclass(pca_result$loadings))
  score$full_score <- as.numeric(pca_result$scores)
  
  #four subcategories scoring
  if(!is.null(set.factor)){
    factor_membership<-set.factor
  }
  for(i in 1:nfactors){
    ex<-subset(x,select=factor_membership[[i]])
    fit<-psych::principal(ex)
    score$ex<-as.numeric(fit$scores)
    names(score)<-c(names(score)[1:(ncol(score)-1)],paste0("factor_",i))
  }
  
  #Step 5: Rescale scores by multiplyin *20, adding 100
  # zip level
  score$full_score <- round(((20*(score$full_score-mean(score$full_score))/sd(score$full_score)) + 100),2)
  nams<-names(score)[(ncol(score)-nfactors+1):ncol(score)]
  for(i in nams){
    score[,i]<-round(((20*(score[,i]-mean(score[,i]))/sd(score[,i])) + 100),2)
  }
  #merge back to geoid labels
  score$rowid<-rownames(score)
  score<-score[order(score$rowid),]
  score<-merge(labs,score,by="rowid",all.x=T)
  
  #Step 6: validate using external data----
  #NOT PERFORMED
  #Step 7: Perform Rankings of overall score at national (or state)----
  score$ranked_score<-rank((score$full_score)/nrow(score),ties.method='min')
  #score$ranked_percentile<-100*((score$full_score-min(score$full_score)))/(max(score$full_score)-min(score$full_score))
  score$ranked_percentile<-100*((score$ranked_score-min(score$ranked_score)))/(max(score$ranked_score)-min(score$ranked_score))
  score$ranked_decile<-ifelse(score$ranked_percentile<10,1,
                              ifelse(score$ranked_percentile<20,2,
                                     ifelse(score$ranked_percentile<30,3,
                                            ifelse(score$ranked_percentile<40,4,
                                                   ifelse(score$ranked_percentile<50,5,
                                                          ifelse(score$ranked_percentile<60,6,
                                                                 ifelse(score$ranked_percentile<70,7,
                                                                        ifelse(score$ranked_percentile<80,8,
                                                                               ifelse(score$ranked_percentile<90,9,
                                                                                      ifelse(score$ranked_percentile>=90,10,NA
                                                                                      ))))))))))
  
  #END----
  output<-list(score=score[,c("NAME",geoid,"full_score","ranked_score","ranked_percentile","ranked_decile")],
               sub_score=score[,c("NAME",geoid,names(score)[(ncol(score)-(2+nfactors)):(ncol(score)-3)])],
               norm.val=x,
               corr.mat=corr_matrix,
               eigen.val=t1,
               eigen.vec=t2,
               scree=p2,
               eigen.cor=var$cos2,
               fa_result=fa_result,
               fa_single=fa_result2,
               factor_members=factor_membership,
               pca_result=pca_result
  )
  output
}

key<-"28bdbbaf6753eed6da2aec6c2eec3157a9f00b39"#will need unique key
tidycensus::census_api_key(key)
ipums.key<-"59cba10d8a5da536fc06b59d66e936288d8a4d7eb286182ea94a5c84"
set_ipums_api_key(api_key=ipums.key)

ADI_construct <- function(comp){
  comp_coef <- as.data.frame(paste0('comp_',1:17))
  comp_coef[,2] <- c(0.0849,-0.0970,-0.0874,-0.0977,0.0936,-0.0688,-0.0781,-0.0770,
                     -0.0615,0.0806,0.0977,0.1037,0.0719,0.0694,0.0877,0.0510,0.0556)
  colnames(comp_coef) <- c('component','value')
  
  temp_score <- as.matrix(comp[,2:18])%*%as.vector(comp_coef[,2])
  scaled_score <- round(((20 * (temp_score-mean(temp_score))/sd(temp_score)) + 100),2)
  
  return(scaled_score)
}

geo_impute<-function(x,geoid="geoid",from,to,type="percent",year=2019){#requires fips coding
  ipums.key<-"59cba10d8a5da536fc06b59d66e936288d8a4d7eb286182ea94a5c84"
  set_ipums_api_key(api_key=ipums.key)
  #hierarchy----
  hierarchy<-factor(c("county","tract","bg"),ordered=T,levels=rev(c("county","tract","bg")))
  #load geographies from ipumsr #tigris----
  geos<-c(from,to)
  geos<-ifelse(geos=="block group","blck_grp",
            ifelse(geos=="zip code tabulation area","zcta",geos))
  files<-get_metadata_nhgis('shapefiles',api_key=ipums.key)
  files<-files[files$year==year&
                 (str_detect(files$name,paste0("us_",geos[1]))|
                    str_detect(files$name,paste0("us_",geos[2]))),]$name
  
  if("county" %in% geos){
    #county<-tigris::counties(year=year)
    extract<-ipumsr::define_extract_nhgis(description=paste0("county_shp_",year),
                                          shapefiles=paste0("us_county_",year,"_tl",year))
    sf<-download_extract(wait_for_extract(submit_extract(extract)))
    county<-read_ipums_sf(sf)
    file.remove(sf)
    remove(sf)
    gc()
    
  }
  if("tract" %in% geos){
    #tract<-NULL
    #for(i in state.abb){
    #  out<-as.data.frame(tigris::tracts(state=i,year=year,resolution="500k"))[,c("STATEFP","GEOID")]
    #  tract<-rbind(tract,out)
    #}
    extract<-ipumsr::define_extract_nhgis(description=paste0("tract_shp_",year),
                                          shapefiles=paste0("us_tract_",year,"_tl",year))
    sf<-download_extract(wait_for_extract(submit_extract(extract)))
    tract<-read_ipums_sf(sf)
    file.remove(sf)
    remove(sf)
    gc()
  }
  if("bg" %in% geos){
    #bg<-NULL
    #for(i in state.abb){
    #  out<-as.data.frame(tigris::block_groups(state=i,year=year))[,c("STATEFP","GEOID")]
    #  bg<-rbind(bg,out)
    #}
    extract<-ipumsr::define_extract_nhgis(description=paste0("blck_grp_shp_",year),
                                          shapefiles=paste0("us_blck_grp_",year,"_tl",year))
    sf<-download_extract(wait_for_extract(submit_extract(extract)))
    bg<-read_ipums_sf(sf)
    file.remove(sf)
    remove(sf)
    gc()
  }
  #create spine on from data----
  if(from=="county"){
    county<-merge(as.data.frame(county[,c("STATEFP","GEOID")]),
                  x,by.x="GEOID",by.y=geoid,all.x=T)
  }
  if(from=="tract"){
    tract<-merge(as.data.frame(tract[,c("STATEFP","GEOID")]),
                 x,by.x="GEOID",by.y=geoid,all.x=T)
  }
  if(from=="bg"){
    bg<-merge(as.data.frame(bg[,"GEOID"]),
              x,by.x="GEOID",by.y=geoid,all.x=T)
  }
  #expand to next level using to geographies-----
  if(from=="county"&to=="tract"){
    tract$tract<-tract$GEOID
    tract$GEOID<-substr(tract$GEOID,1,5)
    county<-merge(county,tract[,c("GEOID","tract")],by="GEOID",all.x=T)#hangs here
    out<-county[,!names(county) %in% c("geometry","STATEFP","GEOID")]
  }
  if(from=="county"&to=="bg"){
    bg$bg<-bg$GEOID
    bg$GEOID<-substr(bg$GEOID,1,5)
    county<-merge(county,bg[,c("GEOID","bg")],by="GEOID",all.x=T)
    out<-county[,!names(county) %in% c("geometry","STATEFP","GEOID")]
  }
  if(from=="tract"&to=="bg"){
    bg$bg<-bg$GEOID
    bg$GEOID<-substr(bg$GEOID,1,11)
    tract<-merge(tract,bg[,c("GEOID","bg")],by="GEOID",all.x=T)
    out<-tract[,!names(tract) %in% c("geometry","STATEFP","GEOID")]
  }
  if(from=="bg"&to=="tract"){
    bg$bg<-bg$GEOID
    bg$GEOID<-substr(bg$GEOID,1,11)
    bg<-merge(bg,tract[,"GEOID"],by="GEOID",all.x=T)
    out<-bg[,!names(bg) %in% c("geometry","STATEFP","GEOID")]
    
  }
  if(from=="bg"&to=="county"){
    bg$bg<-bg$GEOID
    bg$GEOID<-substr(bg$GEOID,1,5)
    bg<-merge(bg,county[,"GEOID"],by="GEOID",all.x=T)
    out<-bg[,!names(bg) %in% c("geometry","STATEFP","GEOID")]
    
  }
  if(from=="tract"&to=="county"){
    tract$tract<-tract$GEOID
    tract$GEOID<-substr(tract$GEOID,1,5)
    tract<-merge(tract,county[,c("GEOID")],by="GEOID",all.x=T)
    out<-tract[,!names(tract) %in% c("geometry","STATEFP","GEOID")]
    
  }
  
  if(type=="percent"|type=="rate"){#if percent going up average, going down impute
    if(hierarchy[hierarchy==from]>hierarchy[hierarchy==to]){#going down
      out<-out[!is.na(out[,to])&!is.na(out[,names(x)[!names(x) %in% geoid]]),]
    }else{#going up
      span<-ifelse(to=="county",5,
                   ifelse(to=="tract",11,
                          ifelse(to=="bg",12)))
      out[,from]<-substr(out[,from],1,span)
      out<-out[!is.na(out[,from])&!is.na(out[,names(x)[!names(x) %in% geoid]]),]
      m<-aggregate(data=out,formula(paste0(names(x)[!names(x) %in% geoid],"~",from)),FUN="mean")
      out<-m
    }
  }
  if(type=="count"){#if count going up sum, going down split
    if(hierarchy[hierarchy==from]>hierarchy[hierarchy==to]){#going down
      span<-ifelse(from=="county",5,
                   ifelse(from=="tract",11,
                          ifelse(from=="bg",12)))
      out$geo<-substr(out[,to],1,span)
      subs<-data.frame(geo=out[,"geo"],val=rep(1,nrow(out)))
      subs<-aggregate(data=subs,val~geo,FUN="sum")
      out<-merge(out,subs,by="geo",all.x=T)
      out<-out[!is.na(out[,to])&!is.na(out[,names(x)[!names(x) %in% geoid]]),]
      out[,names(x)[!names(x) %in% geoid]]<-as.numeric(out[,names(x)[!names(x) %in% geoid]])/out$val
      out<-out[,c(names(x)[!names(x) %in% geoid],to)]
      #      names(out)<-c(names(x)[!names(x) %in% geoid],to)
      out<-out[!duplicated(out),]
      
    }else{#going up
      span<-ifelse(to=="county",5,
                   ifelse(to=="tract",11,
                          ifelse(to=="bg",12)))
      out[,from]<-substr(out[,from],1,span)
      out<-out[!is.na(out[,from])&!is.na(out[,names(x)[!names(x) %in% "geoid"]]),]
      m<-aggregate(data=out,formula(paste0(names(x)[!names(x) %in% "geoid"],"~",from)),FUN="sum")
      out<-m
    }
    
  }
  if(type=="binary"){#if binary going up mode, going down impute binary value
    if(hierarchy[hierarchy==from]>hierarchy[hierarchy==to]){#going down
      out<-out[!is.na(out[,to])&!is.na(out[,names(x)[!names(x) %in% "geoid"]]),]
    }else{#going up
      span<-ifelse(to=="county",5,
                   ifelse(to=="tract",11,
                          ifelse(to=="bg",12)))
      out[,from]<-substr(out[,from],1,span)
      out<-out[!is.na(out[,from])&!is.na(out[,names(x)[!names(x) %in% "geoid"]]),]
      m<-aggregate(data=out,formula(paste0(names(x)[!names(x) %in% "geoid"],"~",from)),FUN="mean")
      m[,names(x)[!names(x) %in% "geoid"]]<-ifelse(m[,names(x)[!names(x) %in% "geoid"]]>0.5,1,0)
      out<-m
    }
    
  }
  
  out
}

score_state<-function(x,geography="county",geoid="GEOID",year=2019){#x must be the score table of makeIndex output
  ipums.key<-"59cba10d8a5da536fc06b59d66e936288d8a4d7eb286182ea94a5c84"
  set_ipums_api_key(api_key=ipums.key)
  geography<-ifelse(geography=="block group","blck_grp",
               ifelse(geography=="zip code tabulation area","zcta",geography))
  files<-get_metadata_nhgis('shapefiles',api_key=ipums.key)
  files<-files[files$year==year&str_detect(files$name,paste0("us_",geography)),]$name
  if(geography=="county"){
    #out<-as.data.frame(tigris::counties(year=year))[,c("STATEFP","GEOID")]
    extract<-ipumsr::define_extract_nhgis(description=paste0("county_shp_",year),
                                          shapefiles=paste0("us_county_",year,"_tl",year))
    sf<-download_extract(wait_for_extract(submit_extract(extract)))
    out<-as.data.frame(read_ipums_sf(sf)[,c("STATEFP","GEOID")])
    file.remove(sf)
    remove(sf)
    gc()
    
  }
  if(geography=="tract"){
    #tract<-NULL
    #for(i in state.abb){
    #  out<-as.data.frame(tigris::tracts(state=i,year=year,resolution="500k"))[,c("STATEFP","GEOID")]
    #  tract<-rbind(tract,out)
    #}
    #out<-tract
    #remove(tract)
    #gc()
    extract<-ipumsr::define_extract_nhgis(description=paste0("tract_shp_",year),
                                          shapefiles=paste0("us_tract_",year,"_tl",year))
    sf<-download_extract(wait_for_extract(submit_extract(extract)))
    out<-as.data.frame(read_ipums_sf(sf)[,c("STATEFP","GEOID")])
    file.remove(sf)
    remove(sf)
    gc()
    
  }
  if(geography=="block group"){
    #bg<-NULL
    #for(i in state.abb){
    #  out<-as.data.frame(tigris::block_groups(state=i,year=year))[,c("STATEFP","GEOID")]
    #  bg<-rbind(bg,out)
    #}
    #out<-bg
    #remove(bg)
    #gc()
    extract<-ipumsr::define_extract_nhgis(description=paste0("blck_grp_shp_",year),
                                          shapefiles=paste0("us_blck_grp_",year,"_tl",year))
    sf<-download_extract(wait_for_extract(submit_extract(extract)))
    out<-as.data.frame(read_ipums_sf(sf)[,c("STATEFP","GEOID")])
    file.remove(sf)
    remove(sf)
    gc()
    
  }
  if(geography=="zip code tabulation area"|geography=="zcta5"|geography=="zcta"){
    #zcta5<-as.data.frame(tigris::zctas(year=year))[,c("STATEFP","GEOID10")]
    zcta5<-as.data.frame(zcta_crosswalk[,c("ZCTA5","GEOID")])
    zcta5$STATEFP<-substr(as.character(zcta5[,"GEOID"]),1,2)
    zcta5[,"GEOID"]<-zcta5[,"ZCTA5"]
    out<-zcta5[,c("STATEFP","GEOID")]
    remove(zcta5)
    gc()
  }
  x<-merge(x,out,by.x=geoid,by.y="GEOID",all.x=T)
  x<-x[!duplicated(x)&!is.na(x$full_score),]
  x$STATEFP<-ifelse(is.na(x$STATEFP),substr(x$GEOID,1,2),x$STATEFP)
  out<-NULL
  for(i in levels(as.factor(x[,"STATEFP"]))){
    o<-x[x[,"STATEFP"]==i,]
    o$state_ranked_score<-o$full_score
    o$state_ranked_percentile<-100*((o$state_ranked_score-min(o$state_ranked_score)))/(max(o$state_ranked_score)-min(o$state_ranked_score))
    o$state_ranked_decile<-ifelse(o$state_ranked_percentile<10,1,
                                ifelse(o$state_ranked_percentile<20,2,
                                       ifelse(o$state_ranked_percentile<30,3,
                                              ifelse(o$state_ranked_percentile<40,4,
                                                     ifelse(o$state_ranked_percentile<50,5,
                                                            ifelse(o$state_ranked_percentile<60,6,
                                                                   ifelse(o$state_ranked_percentile<70,7,
                                                                          ifelse(o$state_ranked_percentile<80,8,
                                                                                 ifelse(o$state_ranked_percentile<90,9,
                                                                                        ifelse(o$state_ranked_percentile>=90,10,NA
                                                                                        ))))))))))
    out<-rbind(out,o)
  }
  out<-merge(x,out[,c(geoid,"state_ranked_percentile","state_ranked_decile")],by=geoid,all.x=T)
  out
}
#score_state(x=obj$score,geography="zcta5",geoid="GEOID",year=2019)

pullACS<-function(geography="county",variables=c("B01001_001"),state=state.abb,geometry=T,year=2019,impute=F,impute.method="rpart",calc="none",calc.adi=F,survey="acs5"){
  ex<-get_acs(geography=geography,variables=variables,state=state,year=year,geometry=geometry,survey=survey)
  newvars<-unique(c(variables,"B01001_001","B25001_001","B26001_001","B19113_001"))
  map<-as.data.frame(ex[,c("GEOID","NAME")])
  ex<-ex[!is.na(ex$estimate),]
  ex<-spread(data=ex[,c("GEOID","variable","estimate")],key=variable,value=estimate,fill=NA)
  #ex[ex==-666666666]<-NA
  #impute base block groups
  if(geography=="block group"&impute==T){
    for(i in newvars[!newvars %in% names(ex)]){
      e<-get_acs(geography="tract",variables=c(i),state=state,year=year,geometry=F,survey=survey)
      e<-spread(data=e[,c("GEOID","variable","estimate")],key=variable,value=estimate,fill=NA)
      if(max(na.omit(e[,i]))==1&min(na.omit(e[,i]))==0){#rate and percent values are hard to interpret, but rely on direct imputation
        type<-"binary"#binary in final form
      }
      if(i %in% c("B01001_001","B25001_001","B26001_001")){#count metrics tend to be skewed and should be split evenly to bg
        type<-"count"#count only, to be split from tract to bg
      }else{#use for central tendency values where direct imputation is required
        type<-"rate"#vary from 1 to 100 in final form
      }
      e<-geo_impute(x=e,geoid="GEOID",from="tract",to="bg",year=year,type=type)
      ex<-merge(ex,e,by.x="GEOID",by.y="bg",all.x=T)
    }
    ex<-ex[,c("GEOID",newvars)]
  }
  base<-ex
  
  if(calc.adi==T){
    adi_plus <- c('B15003_005','B15003_006','B15003_007','B15003_008','B15003_009',
                  'B15003_010','B15003_011','B15003_012','B15002_011','B15002_012',
                  'B15002_013','B15002_014','B15002_015','B15002_016','B15002_017',
                  'B15002_018','B15002_028','B15002_029','B15002_030','B15002_031',
                  'B15002_032','B15002_033','B15002_034','B15002_035','C24050_008',
                  'C24050_009','C24050_010','C24050_011','C24050_012','C24050_013',
                  'C24050_014','B19113_001','B19001_002','B19001_011','B19001_012',
                  'B19001_013','B19001_014','B19001_015','B19001_016','B19001_017',
                  'B25077_001','B25064_001','B25088_002','B25003_002','B23025_005',
                  'B17010_002','C17002_002','B09002_002','B25044_003','B25044_010',
                  'B25043_007','B25043_016','B25016_007','B25016_016','B25014_005',
                  'B25014_006','B25014_007','B25014_011','B25014_012','B25014_013',
                  'B15003_001','C24050_001','B25003_001','B23025_001','B17010_001',
                  'C17002_001','B09002_001','B25044_001','B25043_001','B25016_001','B25014_001',
                  'B08301_001','B08301_010'
                  ) #mod to include prc_publictransit for ICTR projct
    ex<-get_acs(geography=geography,variables=adi_plus,state=state,year=year,geometry=F,survey="acs5")
    ex<-ex[!is.na(ex$estimate),]
    ex<-spread(data=ex[,c("GEOID","variable","estimate")],key=variable,value=estimate,fill=NA)
    #ex[ex==-666666666]<-NA
    raw_dt<-data.table(ex)
    
    if(geography=="block group"){
      select<-c('C24050_001','C24050_008','C24050_009','C24050_010','C24050_011','C24050_012','C24050_013','C24050_014')
      dt_bg_partial <- ex[,!names(ex) %in% select]
      dt_bg_partial$tempid<-substr(dt_bg_partial$GEOID,1,11)
      dt_tract_partial<-get_acs(geography="tract",variables=select,state=state,year=year,geometry=F,survey="acs5")
      dt_tract_partial<-dt_tract_partial[!is.na(dt_tract_partial$estimate),]
      dt_tract_partial<-spread(data=dt_tract_partial[,c("GEOID","variable","estimate")],key=variable,value=estimate,fill=NA)
      outs<-data.table(GEOID=NA)
      for(j in select){
        d<-geo_impute(x=dt_tract_partial[,c("GEOID",j)],geoid="GEOID",from="tract",to="bg",type="count",year=year)
        outs<-merge(outs,d,by.x="GEOID",by.y="bg",all=T)
      }
      outs<-outs[!is.na(outs$GEOID),]
      ### -- impute missing value at bg from tract data
      dt_bg <- merge(dt_bg_partial,outs,by="GEOID",all.x=T)
      dt_bg<-dt_bg[,c("GEOID",adi_plus)]
      raw_dt<-data.table(dt_bg)
      raw_dt<-raw_dt%>%mutate_at(2:ncol(raw_dt),as.numeric)
      raw_dt[is.na(raw_dt) | raw_dt < 0] <- 0
      rm(dt_bg,dt_tract_partial,dt_bg_partial,dt_tract_partial)
      gc()
      
    }
    
    if(geography=="tract"){
      select<-c('C24050_001','C24050_008','C24050_009','C24050_010','C24050_011','C24050_012','C24050_013','C24050_014')
      dt_bg_partial <- ex[,!names(ex) %in% select]
      dt_bg_partial$tempid<-substr(dt_bg_partial$GEOID,1,11)
      dt_tract_partial<-get_acs(geography="tract",variables=select,state=state,year=year,geometry=F,survey="acs5")
      dt_tract_partial<-dt_tract_partial[!is.na(dt_tract_partial$estimate),]
      dt_tract_partial<-spread(data=dt_tract_partial[,c("GEOID","variable","estimate")],key=variable,value=estimate,fill=NA)
      outs<-data.table(GEOID=NA)
      for(j in select){
        d<-geo_impute(x=dt_tract_partial[,c("GEOID",j)],geoid="GEOID",from="tract",to="bg",type="count",year=year)
        outs<-merge(outs,d,by.x="GEOID",by.y="bg",all=T)
      }
      outs<-outs[!is.na(outs$GEOID),]
      ### -- impute missing value at bg from tract data
      dt_bg <- merge(dt_bg_partial,outs,by="GEOID",all.x=T)
      dt_bg<-dt_bg[,c("GEOID",adi_plus)]
      raw_dt<-data.table(dt_bg)
      raw_dt<-raw_dt%>%mutate_at(2:ncol(raw_dt),as.numeric)
      raw_dt[is.na(raw_dt) | raw_dt < 0] <- 0
      rm(dt_bg,dt_tract_partial,dt_bg_partial,dt_tract_partial)
      gc()
      
    }
    
    #Construction of ADI at BG----
    #weights from Singh 2003
    raw_dt<-raw_dt[complete.cases(raw_dt),]
    comp_coef <- as.data.frame(paste0('comp_',1:17))
    comp_coef[,2] <- c(0.0849,-0.0970,-0.0874,-0.0977,0.0936,-0.0688,-0.0781,-0.0770,
                       -0.0615,0.0806,0.0977,0.1037,0.0719,0.0694,0.0877,0.0510,0.0556)
    colnames(comp_coef) <- c('component','value')
    
    #Calculation of scores
    raw_dt[,comp_1:= 100*(B15003_005+B15003_006+B15003_007+B15003_008+B15003_009+B15003_010+B15003_011+B15003_012)/B15003_001]#pctpeoplewithlessthan9thgrade
    raw_dt[,comp_2:= 100*(B15002_011+B15002_012+B15002_013+B15002_014+B15002_015+B15002_016+B15002_017+B15002_018+
                            B15002_028+B15002_029+B15002_030+B15002_031+B15002_032+B15002_033+B15002_034+B15002_035)/B15003_001]#pctpeoplewithatleasthseducation
    raw_dt[,comp_3:= 100*(C24050_008+C24050_009+C24050_010+C24050_011+C24050_012+C24050_013+C24050_014)/C24050_001]#pctwithwhitecollarjob
    raw_dt[,comp_4:= B19113_001]#median family income
    raw_dt[,comp_5:=B19001_002/(B19001_011+B19001_012+B19001_013+B19001_014+B19001_015+B19001_016+B19001_017)]#ratioofthosemakingundert10ktothosemakingover50k
    raw_dt[,comp_5:= 100*log(comp_5)]
    #raw_dt[,comp_5:= 100*comp_5]
    raw_dt[is.infinite(comp_5),comp_5:=0]
    raw_dt[,comp_6:=B25077_001]#medianhousevalue
    raw_dt[,comp_7:=B25064_001]#medianrent
    raw_dt[,comp_8:=B25088_002]#medianmortgage
    raw_dt[,comp_9:=100*B25003_002/B25003_001]#pctowneroccupiedhousing
    raw_dt[,comp_10:=100*B23025_005/B23025_001]#pctpeopleunemployed
    raw_dt[,comp_11:=100*B17010_002/B17010_001]#pctfamiliesinpoverty
    raw_dt[,comp_12:=100*C17002_002/C17002_001]#pctlivingbelow150pctfederalpoverty
    raw_dt[,comp_13:=100*B09002_002/B09002_001]#pcthouseholdswithchildrenwhoaresingleparent
    raw_dt[,comp_14:=100*(B25044_003+B25044_010)/B25044_001]#pcthouseholdswithnovehicles
    raw_dt[,comp_15:=100*(B25043_007+B25043_016)/B25043_001]#pctHouseholdsWithoutTelephone
    raw_dt[,comp_16:=100*(B25016_007+B25016_016)/B25016_001]#pctHouseholdsWithoutPlumbing
    raw_dt[,comp_17:=100*(B25014_005+B25014_006+B25014_007+B25014_011+B25014_012+B25014_013)/B25014_001]#pcthouseholdswithmorethanoneperroom
    
    #Clean output
    comp <- raw_dt[,1]; comp <- cbind(comp,raw_dt[,paste0('comp_',1:17)])
    comp[,2:18][is.na(comp[,2:18])] <- 0
    comp[,2:18][abs(comp[,2:18])>10^10] <- 0
  
    #Create continuous ADI
    comp[,'ADI'] <- ADI_construct(comp)
    comp_out <- comp[,-c(2:18)]
    quantile(comp_out$ADI)
    #plot(density(comp$ADI))
    comp[,'rnk_adi']<-100*percent_rank(comp$ADI)
    #plot(comp$nat_rnk_adi~comp$ADI)#double check: high rank=high deprivation?
    names(comp)<-c("GEOID",
                   "pctpeoplewithlessthan9thgrade",
                   "pctpeoplewithatleasthseducation",
                   "pctwithwhitecollarjob",
                   "medianfamilyincome",
                   "ratioofthosemakingunder10ktothosemakingover50k",
                   "medianhousevalue",
                   "medianrent",
                   "medianmortgage",
                   "pctowneroccupiedhousing",
                   "pctpeopleunemployed",
                   "pctfamiliesinpoverty",
                   "pctlivingbelow150pctfederalpoverty",
                   "pcthouseholdswithchildrenwhoaresingleparent",
                   "pcthouseholdswithnovehicles",
                   "pctHouseholdsWithoutTelephone",
                   "pctHouseholdsWithoutPlumbing",
                   "pcthouseholdswithmorethanoneperroom",
                   "ADI","rnk_adi")
    
    base<-merge(base,comp,by="GEOID",all.x=T)
    
  }
  if(calc=="housing"){
    vars<-data.frame(
      variable=c("housingvalue","mortgage","mortgagecost","mortgagetax","mortgageval","rentgross","occupyown","occupy","rentprcinc","units","unitsownpop","vehicle","fuelheat","kitchen","plumbing","rooms"),
      numerator=c("B25077_001","B25027_002","B25089_002","B25103_001","B25082_002","B25064_001","B25003_002","B25002_002","B25071_001","B25032_001","B25008_002","B25046_001","B25040_002+B25040_003+B25040_004+B25040_005+B25040_008","B25051_002","B25047_002","B25018_001"),
      denominator=c(NA,"B25027_001","B25089_001",NA,"B25082_001",NA,"B25003_001","B25002_001",NA,"B01001_001","B25008_001","B25045_001","B25040_001","B25051_001","B25047_001",NA)
    )
    lvars<-as.character(na.omit(unlist(str_split(unique(c(vars$numerator,vars$denominator)),"[+]"))))
    ex<-get_acs(geography=geography,variables=lvars,state=state,year=year,geometry=F,survey="acs5")
    ex<-ex[!is.na(ex$estimate),]
    ex<-as.data.frame(spread(data=ex[,c("GEOID","variable","estimate")],key=variable,value=estimate,fill=NA))
    ex[ex==-666666666]<-NA
    #remove geographies that are invalid
    ex2<-get_acs(geography=geography,variables=c("B01001_001","B25001_001","B26001_001","B19113_001"),state=state,year=year,geometry=F,survey="acs5")#has population and housing units
    ex2<-tidyr::spread(ex2[,c("GEOID","variable","estimate")],key=variable,value=estimate,fill=0)
    ex2$valid<-ifelse(ex2$B01001_001>100&!is.na(ex2$B01001_001)&#has population >100 per ADI and IMD methods at BG
                        ex2$B25001_001>30&!is.na(ex2$B25001_001)&
                        ex2$B26001_001/ex2$B01001_001<0.333&!is.na(ex2$B26001_001),1,0)#and has housing >30 per ADI and IMD methods at BG
    ex<-merge(ex,ex2[,c("GEOID",names(ex2)[!names(ex2) %in% names(ex)])],by="GEOID",all.x=T)
    ex<-ex[ex$valid==1,]#compute housing on only the valid geographies and removes values for invalid
    remove(ex2)
    gc()
    
    #impute missing raw
    if(impute==T){
      newvars<-unique(c(lvars,"B01001_001","B25001_001","B26001_001","B19113_001"))
      if(geography=="block group"){
        for(i in newvars[!newvars %in% names(ex)]){
          e<-get_acs(geography="tract",variables=c(i),state=state,year=year,geometry=F,survey=survey)
          e<-spread(data=e[,c("GEOID","variable","estimate")],key=variable,value=estimate,fill=NA)
          if(max(na.omit(e[,i]))==1&min(na.omit(e[,i]))==0){#rate and percent values are hard to interpret, but rely on direct imputation
            type<-"binary"#binary in final form
          }
          if(i %in% c("B01001_001","B25001_001","B26001_001")){#count metrics tend to be skewed and should be split evenly to bg
            type<-"count"#count only, to be split from tract to bg
          }else{#use for central tendency values where direct imputation is required
            type<-"rate"#vary from 1 to 100 in final form
          }
          e<-geo_impute(x=e,geoid="GEOID",from="tract",to="bg",year=year,type=type)
          ex<-merge(ex,e,by.x="GEOID",by.y="bg",all.x=T)
        }
        ex<-ex[,c("GEOID",newvars,"valid")]
      }
      
      if(impute.method=="mice"){
        hold<-state
        remove(state)
        gc()
        library(mice)
        ex2<-ex[,colnames(ex)[2:(ncol(ex)-1)]]
        imp<-mice::mice(ex2,maxit=2,m=5,seed=1,method="cart")#does not run if there is anything in the environment called "state"
        ex2<-mice::complete(imp,action="long")
        ex2<-ex2[,2:ncol(ex2)]%>%
          group_by(.id)%>%
          summarise_each(funs=c("mean"))
        ex2<-as.data.frame(ex2)
        #remove(ex2)
        gc()
        #imputing one value from ten iterations using a cart to handle univariate values
        state<-hold
        remove(hold)
        gc()
        
        ex$rowid<-rownames(ex)
        #impute<-NULL
        #for(i in names(imp[["imp"]])){
        #  if(length(imp[["imp"]][[i]][[1]])>0){
        #    impute<-c(impute,i)
        #  }
        #}
        impute<-NULL
        for(i in names(ex)){
          if(length(summary(as.factor(is.na(ex[,i]))))>1){
            impute<-c(impute,i)
          }
        }
        
        for(i in impute){
          #e<-data.frame(rowid=rownames(imp[["imp"]][[i]]),X1=imp[["imp"]][[i]][,1])
          e<-ex2[,c(".id",i)]
          names(e)<-c("rowid","X1")
          ex<-merge(ex,e,by="rowid",all.x=T)
          ex[,i]<-ifelse(is.na(ex[,i])&ex[,"valid"]==1,
                         ifelse(is.na(ex$X1),mean(na.omit(ex$X1)),ex$X1),#impute mice or mean if mice is still an NA value
                         ex[,i])
          ex<-ex[,!names(ex) %in% c("X1")]
        }
        ex<-ex[,!names(ex) %in% c("rowid","valid")]

      }else{
        miss<-NULL
        imputed<-NULL
        for(i in colnames(ex)[2:(ncol(ex)-1)]){#model imputation of raw attributes
          if(impute.method=="rpart"){
            f<-paste0(i,"~",paste(c(newvars[!newvars %in% i]),collapse="+"))#at minimum, will predict using pop,housing and group housing stats
            #e<-ex[,c("GEOID",c(newvars[!newvars %in% i]))]
            e<-ex[,c("GEOID",newvars)]
            for(j in names(e)){#mean imputation on our training data
              e[,j]<-ifelse(is.na(e[,j]),mean(na.omit(e[,j])),e[,j])#impute mean for each individual attribute missing values in e
            }
            m<-rpart::rpart(formula=f,data=e,method="anova")#use partition to return values of each attribute using e to fit
            expected<-predict(m,newdata=ex)#store expected values from raw data
            print(paste0(i," model variance explained ",
                         round(1-(sum((expected[!is.na(ex[,i])]-ex[!is.na(ex[,i]),i])^2)/sum(ex[!is.na(ex[,i]),i]^2)),3)))
            #e[,i]<-ifelse(ex$valid==0,NA,#if not valid, return NA
            #               ifelse(is.na(e[,i]),expected,e[,i]))#if variable is na return imputed value
            if(1-(sum((expected[!is.na(ex[,i])]-ex[!is.na(ex[,i]),i])^2)/sum(ex[!is.na(ex[,i]),i]^2))<0.2&
               !is.na(1-(sum((expected[!is.na(ex[,i])]-ex[!is.na(ex[,i]),i])^2)/sum(ex[!is.na(ex[,i]),i]^2)))){#r2 is na if expected observed is one value only (all)
              miss<-c(miss,i)
            }
            imputed<-cbind(imputed,expected)
          }
          if(impute.method=="lm"){
            f<-paste0(i,"~",paste(c(newvars[!newvars %in% i]),collapse="+"))#at minimum, will predict using pop,housing and group housing stats
            #e<-ex[,c("GEOID",c(newvars[!newvars %in% i]))]
            e<-ex[,c("GEOID",newvars)]
            for(j in names(e)){#mean imputation on our training data
              e[,j]<-ifelse(is.na(e[,j]),mean(na.omit(e[,j])),e[,j])#impute mean for each individual attribute missing values in e
            }
            m<-lm(e,formula=f)
            expected<-predict(m,newdata=e)#store expected values from mean-imputed data
            print(paste0(i," model variance explained ",
                         round(1-(sum((expected[!is.na(ex[,i])]-ex[!is.na(ex[,i]),i])^2)/sum(ex[!is.na(ex[,i]),i]^2)),3)))
            if(1-(sum((expected[!is.na(ex[,i])]-ex[!is.na(ex[,i]),i])^2)/sum(ex[!is.na(ex[,i]),i]^2))<0.2&
               !is.na(1-(sum((expected[!is.na(ex[,i])]-ex[!is.na(ex[,i]),i])^2)/sum(ex[!is.na(ex[,i]),i]^2)))){#r2 is na if expected observed is one value only (all)
              miss<-c(miss,i)
            }
            imputed<-cbind(imputed,expected)
            
          }
        }
        colnames(imputed)<-colnames(ex)[2:(ncol(ex)-1)]
        for(i in colnames(ex)[2:(ncol(ex)-1)]){
          if(!i %in% c(miss,"B01001_001","B25001_001","B26001_001","B19113_001")){#skip any attributes where r2 is very poor (<.2) or one of my starting predictors
            gmin<-min(na.omit(ex[,i]))
            gmax<-max(na.omit(ex[,i]))
            ex[,i]<-ifelse(is.na(ex[,i]),
                           imputed[,i],
                           #e[,i],
                           ex[,i])
            ex[,i]<-ifelse(ex[,i]<gmin,gmin,ifelse(ex[,i]>gmax,gmax,ex[,i]))
          }
        }
        ex<-ex[,c("GEOID",lvars)]
      }
    }
    
    for(i in vars$variable){
      #print(i)
      num<-0
      if(stringr::str_count(vars[vars$variable==i,]$numerator,"[+]")>0){
        for(j in 1:(stringr::str_count(vars[vars$variable==i,]$numerator,"[+]")+1)){
          #id<-vars[vars$variable==i,]$variable
          id<-stringr::str_split(vars[vars$variable==i,]$numerator,"[+]")[[1]][j]
          num<-num+as.numeric(ex[,id])
        }
      }else{
        num<-as.numeric(ex[,vars[vars$variable==i,]$numerator])
      }
      
      den<-0
      if(!is.na(vars[vars$variable==i,]$denominator)){
        if(stringr::str_count(vars[vars$variable==i,]$denominator,"[+]")>0){
          for(j in 1:(stringr::str_count(vars[vars$variable==i,]$denominator,"[+]")+1)){
            #id<-vars[vars$variable==i,]$variable
            id<-stringr::str_split(vars[vars$variable==i,]$denominator,"[+]")[[1]][j]
            den<-den+as.numeric(ex[,id])
          }
        }
        else{
          den<-as.numeric(ex[,vars[vars$variable==i,]$denominator])
          den[is.na(den)]<-0
        }
      }
      
      if(length(den)==1){#THIS LINE WAS EDITED TO CORRECT MISMATCH EVALUATION 2-4-25
        ex[,i]<-num
      }else{
        ex[,i]<-num/den
      }
    }
    #clean
    ex<-ex[,c("GEOID",vars$variable)]
    
    ex$fuelheat<-100*ex$fuelheat
    ex$kitchen<-ifelse(ex$kitchen>1&!is.na(ex$kitchen),100,100*ex$kitchen)
    ex$mortgage<-100*ex$mortgage
    ex$mortgageval<-ifelse(ex$mortgageval>1&!is.na(ex$mortgageval),100,100*ex$mortgageval)
    ex$occupy<-100*ex$occupy
    ex$occupyown<-ifelse(ex$occupyown>1&!is.na(ex$occupyown),100,100*ex$occupyown)
    ex$plumbing<-100*ex$plumbing
    ex$units<-100*ex$units
    ex$unitsownpop<-100*ex$unitsownpop
    #ex$phone<-100*ex$phone
    ex$mortgagecost<-ifelse(ex$mortgagecost>1&!is.na(ex$mortgagecost),100,100*ex$mortgagecost)
    ex$vehicle<-ifelse(is.infinite(ex$vehicle),NA,ifelse(ex$vehicle>5,5,ex$vehicle))
    ex$NAME<-as.character(1:nrow(ex))
    out<-makeIndex(x=ex[complete.cases(ex),],geoid="GEOID",nfactors=4,state=NULL,
                   set.factor=list(financial=c("housingvalue","mortgage","mortgagecost","mortgagetax","mortgageval","rentgross"),
                                   ownership=c("occupyown","unitsownpop","vehicle","rooms"),
                                   building_quality=c("fuelheat","kitchen","plumbing"),
                                   other=c("occupy","rentprcinc","units")
                   ))
    housing<-merge(ex[,c("GEOID",vars$variable)],
                   out$score[,c("GEOID","full_score","ranked_percentile")],
                   by="GEOID",all.x=T)
    housing<-merge(housing,
                   out$sub_score[,c("GEOID","factor_1","factor_2","factor_3","factor_4")],
                   by="GEOID",all.x=T)
    housing<-score_state(x=housing,geography=geography,geoid="GEOID",year=year)
    housing<-housing[,!names(housing) %in% c("STATEFP","state_ranked_decile")]
    names(housing)<-c("GEOID",vars$variable,"housing_score","housing_rnk","financial_score","ownership_score","building_quality_score","other_housing_score",
                      "housing_state_rnk")
    
    base<-merge(base,housing,by="GEOID",all.x=T)
    gc()
  }
  base<-merge(base,map[,!names(map) %in% c("geometry")],by="GEOID",all.x=T)
  base[!duplicated(base),]
}

testing<-F
if(testing==T){
  head(pullACS(geography="county",variables=c("B01001_001","B01002_001"),state=state.abb,geometry=T,year=2019,calc.adi=T,survey="acs5",calc="housing"))
  
  head(pullACS(geography="zip code tabulation area",variables=c("B01001_001","B01002_001"),state=state.abb,geometry=T,year=2019,calc.adi=T,survey="acs5",calc="housing"))

  exs<-NULL
  for(j in c("county","zip code tabulation area","tract")){
    saveRDS(pullACS(geography=j,state=state.abb,geometry=T,year=2019,calc="housing",calc.adi=F,survey="acs5"),
            paste0("C:/Users/chris/OneDrive/Desktop/GeoHealth/data/housing_",j,"_2019_sf.rds"))
    
    ex<-pullACS(geography=j,state=state.abb,year=2019,calc="housing",calc.adi=F,geometry=F,survey="acs5")
    ex<-pullACS(geography=j,state=state.abb,year=2019,calc="housing",calc.adi=F,geometry=F,survey="acs5",impute=T)
    
    nams<-names(ex)[!names(ex) %in% c("B01001_001","GEOID","housing_score","housing_rnk","financial_score","ownership_score","building_quality_score","other_housing_score")]
    for(i in nams){
      nam<-paste0("z_",i)
      #ex[,nam]<- liver::transform(ex[,i], method = "zscore") #x is now normalized THIS STOPPED WORKING
      ex[,nam]<-liver::zscore(ex[,i])
    }
    
    ex$NAME<-as.character(1:nrow(ex))
    ex$geography<-j
    exs<-rbind(exs,ex)
    obj<-makeIndex(x=ex[,names(ex) %in% c("GEOID","NAME",nams)],geoid="GEOID",nfactors=4,state=NULL,
                   set.factor=list(financial=c("housingvalue","mortgage","mortgagecost","mortgagetax","mortgageval","rentgross"),
                                   ownership=c("occupyown","unitsownpop","vehicle","rooms"),
                                   building_quality=c("fuelheat","kitchen","plumbing"),
                                   other=c("occupy","rentprcinc","units"))
    )
    saveRDS(obj,paste0("C:/Users/chris/OneDrive/Desktop/GeoHealth/data/housing_",j,"_2019_output.rds"))
    
  }
  saveRDS(exs[,c("GEOID","geography",nams,paste0("z_",nams),"housing_score","housing_rnk","financial_score","ownership_score","building_quality_score","other_housing_score")],
          paste0("C:/Users/chris/OneDrive/Desktop/GeoHealth/data/housing_2019.rds"))
}

if(testing==T){
  test1<-pullACS(geography="zip code tabulation area",state=state.abb,geometry=F,year=2019,impute=F,calc="housing")
  test2<-pullACS(geography="zip code tabulation area",state=state.abb,geometry=F,year=2019,impute=T,impute.method="rpart",calc="housing")
  #test2<-pullACS(geography="zip code tabulation area",state=state.abb,geometry=F,year=2019,impute=T,impute.method="lm",calc="housing")
  #test2<-pullACS(geography="zip code tabulation area",state=state.abb,geometry=F,year=2019,impute=T,impute.method="mice",calc="housing")#mice broken inside of pullACS function
  source("C:/Users/chris/OneDrive/Desktop/GeoHealth/scripts/jaccard.r")
  
  jaccard(round(test1[,2:(ncol(test1)-1)],1),round(test2[rownames(test1),2:(ncol(test2)-1)],1))
  summary(test1)
  summary(test2)
  summary(test1[,2:(ncol(test1)-1)]-test2[rownames(test1),2:(ncol(test2)-1)])
  
}
