#retrieve our shp files from tigris----
library(sf)#open shp files
library(tigris)#identifies zip by state and returns shp!!!
#load tmap and leaflet----
library(tmap)#main plotting library
library(leaflet)#interactive plotting for multiple layers and aesthetics
library(RColorBrewer)#adds color palettes to environment

flat_map<-function(data=NULL,#should be a data.frame with a valid geoid for fips or zcta5
                   shp=NULL,#can be manually provided but must be type sf data.frame from .shp (see Census)
                   year=2010,#should coincide with decennial Census and reflects updates to geometries
                   state=NULL,#most relevant for plotting zip, tract and bg. Please specify for these
                   geography="county",#select geography: state, county, zcta5, tract or bg
                   geoid=NULL,#user specified as whatever variable name exists in data file
                   var=NULL,#can be a single value for flat plotting or vector of variable names in data
                   type="flat",#specifies whether to return flat map, leaflet or multiple layer leaflet
                   bin=5,#number of levels to each fill in plot
                   palette="RdBu"#RColorBrewer palette for plotting fills. See brewer.pal.info
){
  if(is.null(shp)){#if shp is not specified, retrieve defaults from tigris
    if(geography=="state"){
      shp<-states(year=year)#ignores state argument
    }
    if(geography=="county"){#select national counties or state specific
      if(is.null(state)){
        shp<-counties(year=year)
      }else{
        shp<-counties(state=state,year=year)
      }
    }
    #no choice on smaller geographies, must have a state argument
    if(geography=="zcta5"){
      shp<-zctas(state=state,year=year)#wait this is shp and works so much better!!!
    }
    if(geography=="tract"){
      shp<-tracts(state=state,year=year)#ok tigris can do it all?
    }
    if(geography=="bg"){
      shp<-block_groups(state=state,year=year)
    }
  }
  
  #the function merges and plots values specified with bin and palette being the only adjusted aesthetics
  if("GEOID10" %in% names(as.data.frame(shp))){#if not geoid10, then error
    if("ZCTA5CE10" %in% names(as.data.frame(shp))){#if zip, then use zcta5ce10
      shp<-merge(shp,data,by.x="ZCTA5CE10",by.y=geoid,all.x=T)
    }else{#if not zip, then use geoid10
      shp<-merge(shp,data,by.x="GEOID10",by.y=geoid,all.x=T)
    }
  }else{#return error, no shp
    print("Error: No suitable geoid on shp file. Must use Census provided geometries on tigris package or similar .shp file")
  }
  if(length(var)==1){#for single variable plotting
    tm<-tm_shape(shp)+
      tm_fill(var,n=bin,
              palette=brewer.pal(bin,palette),
              textNA="No Data",
              colorNA="black",
              alpha=0.55,
              border.col="gray20",
              border.alpha=0.8,
              group=var,
              midpoint=NA)+
      tm_borders()+
      tm_scale_bar()
  }else{
    tm<-NULL
    for(i in 1:length(var)){
      tmx<-tm_shape(shp)+
        tm_fill(var[i],n=bin,
                palette=brewer.pal(bin,palette),
                textNA="No Data",
                colorNA="black",
                alpha=0.55,
                border.col="gray20",
                border.alpha=0.8,
                group=var[i],
                midpoint=NA)+
        tm_borders()+
        tm_scale_bar()
      tm<-tm+tmx
    }
  }
  if(type=="leaflet"|type=="multi"){
    tmap_leaflet(tm,show=F)
  }else{
    tm
  }
  #  tm
}

#the loadSHP function should streamline this a bit----
loadSHP<-function(state=NULL,geography="tract",year="2020",outdir=NULL,destfile=NULL){
  st<-tigris::fips_codes
  st<-st[st$state==state,]$state_code[[1]]#converts state abbreviation to fips code
  
  if(is.null(state)){#if whole nation, only retrieve .shp for state/county
    if(geography=="state"){#for state
      url<-paste0("https://www2.census.gov/geo/tiger/TIGER",year,"/STATE/tl_",year,"_us_state.zip")
    }
    if(geography=="county"){
      url<-paste0("https://www2.census.gov/geo/tiger/TIGER",year,"/COUNTY/tl_",year,"_us_county.zip")
    }
    if(geography=="zcta5"){
      url<-paste0("https://www2.census.gov/geo/tiger/TIGER",year,"/ZCTA510/tl_",year,"_us_zcta510.zip")
    }
    
  }else{#if state was specified, retrieve only for that area
    if(geography=="county"){#still pulls everything if county, which should not be an issue
      url<-paste0("https://www2.census.gov/geo/tiger/TIGER",year,"/COUNTY/tl_",year,"_us_county.zip")
    }
    if(geography=="zcta5"){
      url<-paste0("https://www2.census.gov/geo/tiger/TIGER",year,"/ZCTA510/tl_",year,"_us_zcta510.zip")
    }
    if(geography=="tract"){#still pulls everything if zip, which should not be an issue
      url<-paste0("https://www2.census.gov/geo/tiger/TIGER",year,"/TRACT/tl_",year,"_",st,"_tract.zip")
    }
    if(geography=="bg"){
      url<-paste0("https://www2.census.gov/geo/tiger/TIGER",year,"/BG/tl_",year,"_",st,"_bg.zip")
    }
  }
  
  download.file(url=url,destfile=destfile)#this retrieves from Census
  unzip(zipfile=destfile,exdir=outdir)#unzip to new output directory
  shp_files<-dir(outdir)[substr(dir(outdir),nchar(dir(outdir))-3,nchar(dir(outdir)))==".shp"]
  shp_files<-shp_files[str_detect(url,substr(shp_files,1,nchar(shp_files)-4))]
  
  shp<-st_read(dsn=paste0(outdir,paste0("\\",shp_files)))#then load directly as shp assuming only 1 file matches
  shp
}
