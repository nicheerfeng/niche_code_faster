

###############################################################
########## 暂时保存：##########################################
###############################################################

## ppp_spatstat--点转ppp模型--空间生态学：####
ppp_spatstat <- function(x){
  occ <- x
  cord.dec = SpatialPoints(occ, 
                           proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  )
  occ1 <- spTransform(cord.dec, CRS("+init=epsg:3857"))
  pts <- occ1@coords
  mcp <- function (xy) {
    library(rgdal)
    xy <- as.data.frame(coordinates(xy))
    coords.t <- chull(xy[, 1], xy[, 2])
    xy.bord <- xy[coords.t, ]
    xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
  }
  ranges <- mcp(occ) %>% rgeos::gBuffer(width =2)
  crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  
  proj4string(ranges) <- crs.geo  
  rangee <- spTransform(ranges, CRS("+init=epsg:3857"))
  library(maptools) 
  rangeowin <- as.owin(rangee)
  p <- ppp(pts[,1], pts[,2], window=rangeowin )
  return(p)
}

## zero_1--数据0-1标准化；####
zero_1 <- function(x){
  x <- data.frame(x)
  col <- dim(x)[2]
  for(i in 1:col){
    x[,col+i] <-  (x[,i]-min(x[,i]))/(max(x[,i])-min(x[,i]))
  }
  return(x)
}

## ma_split--机器学习建模数据快速分割；####
ma_split <- function(data,SplitRatio2 = 0.7){
  library(caTools)
  data2 <- cbind(1:dim(data)[1],data) %>% data.frame()
  names(data2)[1] <- "namexx"
  split <- sample.split(data2$namexx, SplitRatio = 0.7)
  # split data set
  train.set <- data[split==TRUE,]
  val.set <- data[split==FALSE,]
  out <- list(train.set,val.set)
  names(out) <- c("train.set","val.set")
  return(out)
}
out <- ma_split(scau)







###############################################################
############# 第一部分： 物种分布数据清理 #####################
###############################################################
## thin4——分布数据去重、去缺失，空间重抽样 一条龙服务 ####
thin4 <- function(x,y){
  occ <- x
  names(occ) <- c("longitude","latitude")
  occs.dups <- duplicated(occ)
  occs <- occ[!occs.dups,]
  # remove NAs
  occs <- occs[complete.cases(occs$longitude, occs$latitude), ]
  occs <- cbind(c(rep("occ",dim(occs)[1])),occs)
  names(occs)[1] <- "name"
  # 空间清理：
  thin=y
  output <- spThin::thin(occs, 
                         'latitude', 'longitude', 'name',
                         thin.par = thin, reps = 100, locs.thinned.list.return = TRUE, 
                         write.files = FALSE, verbose = FALSE)
  # find the iteration that returns the max number of occurrences
  maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
  # if there's more than one max, pick the first one
  maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
  # subset occs to match only thinned occs
  occs <- occs[as.numeric(rownames(maxThin)),]
  return(occs)
  
}

## cleandata--物种分布数据快速清理 ####
cleanocc <- function(x,y){ ## x为经纬度数据，y为筛选的分辨率；
  occn <- x
  names(occn) <- c('longitude', 'latitude')
  occs.dups <- duplicated(occn[c('longitude', 'latitude')])
  occs <- occn[!occs.dups,]
  occs <- occs[complete.cases(occs$longitude, occs$latitude), ]
  occs <- rep("occ",dim(occs)[1]) %>% cbind(.,occs) %>% data.frame()
  names(occs)[1] <- "name"
  output <- spThin::thin(occs, 
                         'latitude', 'longitude', 'name', 
                         thin.par = y, reps = 1, locs.thinned.list.return = TRUE, 
                         write.files = FALSE, verbose = FALSE)
  # find the iteration that returns the max number of occurrences
  maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
  # if there's more than one max, pick the first one
  maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
  # subset occs to match only thinned occs
  occs <- occs[as.numeric(rownames(maxThin)),]  
  return(occs)
}





## variogram2--检验空间自相关- 空间半方差估计 #####
variogram2 <- function(occ,env){ ##x 为分布经纬度；y对应的环境信息；
  occ189vx <- raster::extract(env,occ) %>%  cbind(occ) %>% na.omit() %>% data.frame()
  occ189pr <- point_wgs_utm(occ)
  
  endtry <- cbind(occ189pr@coords,occ189vx)
  coordinates(endtry) = ~longitude +latitude
  
  c.variog <- variogram(deparse(substitute(env)) ~1, endtry, xlab= "Distance (m)", ylab="Semivariance")
  plot(c.variog,main="4")

}








## point_wgs_utm--分布点由地理坐标系转为投影坐标系：####
point_wgs_utm <- function(x){
  occ <- x
  cord.dec = SpatialPoints(occ, 
                           proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  )
  occ <- spTransform(cord.dec, CRS("+init=epsg:3857"))
  return(occ)
}

## poly_mcp_wgs_utm--点构建最小面积并转成投影坐标系:####
poly_mcp_wgs_utm <- function(x){
  mcp <- function (xy) {
    library(rgdal)
    xy <- as.data.frame(coordinates(xy))
    coords.t <- chull(xy[, 1], xy[, 2])
    xy.bord <- xy[coords.t, ]
    xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
  }
  occ <- x
  ranges <- mcp(occ) %>% rgeos::gBuffer(width =2)
  crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  
  proj4string(ranges) <- crs.geo  
  rangee <- spTransform(ranges, CRS("+init=epsg:3857"))
  return(rangee)
}






###############################################################
########### 第二部分：环境数据清理 ############################
###############################################################
## mcp--构建范围 mcp:####
mcp <- function (xy) {
  library(rgdal)
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
}

## frange--构建范围 bound：####
frange <- function(x){
  occs <- x
  xmin <- min(occs$longitude)
  xmax <- max(occs$longitude)
  ymin <- min(occs$latitude)
  ymax <- max(occs$latitude)
  bb <- matrix(c(xmin-4, xmin-4, xmax+4, xmax+4, xmin-4, ymin-3, ymax+4, ymax+4, ymin-3, ymin-3), ncol=2)
  bgExt <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(bb)), 1)))
  return(bgExt)
}

## buffer0_25--构建范围缓冲区 ####
buffer0_25 <- function(x,n){ ## x为经纬度数据，n为缓冲区的精度；
  occs <- x
  occs.xy <- occs
  sp::coordinates(occs.xy) <- ~ longitude + latitude
  bgExt <- occs.xy
  bgExt3 <- rgeos::gBuffer(bgExt,width = n)
  library(maptools)
  data("wrld_simpl")
  s1 <- raster::intersect(wrld_simpl,bgExt3)
  bg_xy2 <- sp::spsample(s1, n, type="regular")
  return(as.data.frame(bg_xy2@coords))
}


## kde_bg--基于物种分布密度构建缓冲区，构建核密度分布：####
kde_bg <- function(x,n){
  library(adehabitatHR)
  as2 = SpatialPoints(cbind(x$longitude, x$latitude), proj4string = CRS("+init=epsg:4326" ))
  as.sp <- spTransform(as2, CRS("+init=epsg:3857"))
  kernel.ref <- kernelUD(as.sp, h = "href", grid = 1000,hlim = c(0.5, 1.5))  # href = the reference bandwidth
  as.kernel.poly <- getverticeshr(kernel.ref, percent = 90) 
  zk_area<- spTransform(as.kernel.poly, CRS("+init=epsg:4326"))
  library(maptools)
  data("wrld_simpl")
  s1 <- raster::intersect(wrld_simpl,zk_area)
  bg_xy2 <- sp::spsample(s1, n, type="regular")
  return(as.data.frame(bg_xy2@coords))
}

kde_range <- function(x){
  library(adehabitatHR)
  as2 = SpatialPoints(cbind(x$longitude, x$latitude), proj4string = CRS("+init=epsg:4326" ))
  as.sp <- spTransform(as2, CRS("+init=epsg:3857"))
  kernel.ref <- kernelUD(as.sp, h = "href", grid = 1000,hlim = c(0.5, 1.5))  # href = the reference bandwidth
  as.kernel.poly <- getverticeshr(kernel.ref, percent = 90) 
  zk_area<- spTransform(as.kernel.poly, CRS("+init=epsg:4326"))
  return(zk_area)
}


## outshp--将空间数据shp转为空间数据框：SpatialPolygons  to SpatialPolygonsDataFrame ####
# 转成可导出的形式：
outshp <- function(x){
  bc <- x
  df <- data.frame(ID=character(), stringsAsFactors=FALSE )
  for (i in bc@polygons ) { df <- rbind(df, data.frame(ID=i@ID, stringsAsFactors=FALSE))  }
  # and set rowname=ID
  row.names(df) <- df$ID
  ##############################################
  rangepy <- SpatialPolygonsDataFrame(bc, df)
  plot(rangepy)
  return(rangepy)
}

## raster_utm--栅格转投影：####
raster_utm <- function(x,y){
  mcp <- function (xy) {
    library(rgdal)
    xy <- as.data.frame(coordinates(xy))
    coords.t <- chull(xy[, 1], xy[, 2])
    xy.bord <- xy[coords.t, ]
    xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
  }
  occ <- x
  env <- stack(y)
  ranges <- mcp(occ) %>% rgeos::gBuffer(width =2)
  env <- mask(crop(env,ranges),ranges)
  crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  
  projection(env) <- crs.geo
  crs2 <-  CRS('+init=EPSG:3857')
  pz <-  projectRaster(env, crs=crs2)
  return(pz)
}


### ned_ret--根据研究区的分布面积确定随机点数量 ####
ned_ret <- function(xy){
  occ <- xy
  
  mcp <- function (xy) {
    library(rgdal)
    xy <- as.data.frame(coordinates(xy))
    coords.t <- chull(xy[, 1], xy[, 2])
    xy.bord <- xy[coords.t, ]
    xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
  }
  
  rang1 <- mcp(occ) 
  rang2 <- rgeos::gBuffer(bgExt,width = 0.1)
  area1 <- raster::area(rang2) * 0.000001
  ned <- function(yy){
    needle = 0
    x = yy
    if(x<100){
      needle = 1000
      
    }else if(x>100 && x<500){
      needle = 2000
      
    }else if(x>500 && x<1000){
      needle = 3000
      
    }else if(x>1000 && x<2000){
      needle = 4000
    }else{needle = 5000
    }
    return(needle)
  }
  ned_va = ned(area1)
  return(ned_va)
  
}






###############################################################
######### 第三部分：可视化快捷命令 ####
###############################################################

## mapview2-- 分布点地图可视化：####
mapview2 <- function(occ){ ## 注意需要保留数据的经纬度信息为：longitude+latitude形式；
  library(mapview)
  # install.packages("mapview")
  library(sp)
  coordinates(occ) <- ~longitude+latitude
  proj4string(occ) <- CRS("+init=epsg:4326")
  # plot data
  mapView(occ)
}
## plotlm--二维变量快速线性相关性统计可视化####
plotlm <- function(x,y){
  xlab2 <- deparse(substitute(x))
  ylab2 <- deparse(substitute(y))
  library(basicTrendline)
  trendline(x,y, model="line2P", ePos.x = "topleft", summary=TRUE,
            eDigit=5,xlab = xlab2,ylab= ylab2)
}
# plotlm(iris$Sepal.Length,iris$Sepal.Width)
## marginal_plot--二维分组数据可视化，将数据集分组后，添加线性统计及边缘密度图；####
## 参考自：https://github.com/ChrKoenig/R_marginal_plot
marginal_plot = function(x, y, group = NULL, data = NULL, lm_show = FALSE, lm_formula = y ~ x, bw = "nrd0", adjust = 1, alpha = 1, plot_legend = T, ...){
  require(scales)
  ###############
  # Plots a scatterplot with marginal probability density functions for x and y. 
  # Data may be grouped or ungrouped. 
  # For each group, a linear fit can be plotted. It is hidden by default, but can be shown by providing lm_show = TRUE.
  # The model can be modified using the 'lm_formula' argument.
  # The 'bw' and 'adjust' argument specify the granularity used for estimating probability density functions. See ?density for more information.
  # For large datasets, opacity may be decreased by setting alpha to a value between 0 and 1.
  # Additional graphical parameters are passed to the main plot, so you can customize axis labels, titles etc.
  ###############
  moreargs = eval(substitute(list(...)))
  
  # prepare consistent df
  if(missing(group)){
    if(missing(data)){
      if(length(x) != length(y)){stop("Length of arguments not equal")}
      data = data.frame(x = as.numeric(x), y = as.numeric(y))
    } else {
      data = data.frame(x = as.numeric(data[,deparse(substitute(x))]), 
                        y = as.numeric(data[,deparse(substitute(y))]))
    }
    if(sum(!complete.cases(data)) > 0){
      warning(sprintf("Removed %i rows with missing data", sum(!complete.cases(data))))
      data = data[complete.cases(data),]
    }
    group_colors = "black"
  } else {
    if(missing(data)){
      if(length(x) != length(y) | length(x) != length(group)){stop("Length of arguments not equal")}
      data = data.frame(x = as.numeric(x), y = as.numeric(y), group = as.factor(group))
    } else {
      data = data.frame(x = as.numeric(data[,deparse(substitute(x))]), 
                        y = as.numeric(data[,deparse(substitute(y))]),
                        group = as.factor(data[,deparse(substitute(group))]))
    }
    if(sum(!complete.cases(data)) > 0){
      warning(sprintf("Removed %i rows with missing data", sum(!complete.cases(data))))
      data = data[complete.cases(data),]
    }
    data = subset(data, group %in% names(which(table(data$group) > 5)))
    data$group = droplevels(data$group)
    group_colors = rainbow(length(unique(data$group)))
  } 
  
  # log-transform data (this is need for correct plotting of density functions)
  if(!is.null(moreargs$log)){
    if(!moreargs$log %in% c("y", "x", "yx", "xy")){
      warning("Ignoring invalid 'log' argument. Use 'y', 'x', 'yx' or 'xy.")
    } else {
      data = data[apply(data[unlist(strsplit(moreargs$log, ""))], 1, function(x) !any(x <= 0)), ]
      data[,unlist(strsplit(moreargs$log, ""))] = log10(data[,unlist(strsplit(moreargs$log, ""))])
    }
    moreargs$log = NULL # remove to prevent double logarithm when plotting
  }
  
  # Catch unwanted user inputs
  if(!is.null(moreargs$col)){moreargs$col = NULL}
  if(!is.null(moreargs$type)){moreargs$type = "p"}
  
  # get some default plotting arguments
  if(is.null(moreargs$xlim)){moreargs$xlim = range(data$x)} 
  if(is.null(moreargs$ylim)){moreargs$ylim = range(data$y)}
  if(is.null(moreargs$xlab)){moreargs$xlab = deparse(substitute(x))}
  if(is.null(moreargs$ylab)){moreargs$ylab = deparse(substitute(y))}
  if(is.null(moreargs$las)){moreargs$las = 1} 
  
  # plotting
  tryCatch(expr = {
    ifelse(!is.null(data$group), data_split <- split(data, data$group), data_split <- list(data))
    orig_par = par(no.readonly = T)
    par(mar = c(0.25,5,1,0))
    layout(matrix(1:4, nrow = 2, byrow = T), widths = c(10,3), heights = c(3,10))
    
    # upper density plot
    plot(NULL, type = "n", xlim = moreargs$xlim, ylab = "density",
         ylim = c(0, max(sapply(data_split, function(group_set) max(density(group_set$x, bw = bw)$y)))), main = NA, axes = F)
    axis(2, las = 1)
    mapply(function(group_set, group_color){lines(density(group_set$x, bw = bw, adjust = adjust), col = group_color, lwd = 2)}, data_split, group_colors)
    
    # legend
    par(mar = c(0.25,0.25,0,0))
    plot.new()
    if(!missing(group) & plot_legend){
      legend("center", levels(data$group), fill = group_colors, border = group_colors, bty = "n", title = deparse(substitute(group)), title.adj = 0.1)
    }
    
    # main plot
    par(mar = c(4,5,0,0))
    if(missing(group)){
      do.call(plot, c(list(x = quote(data$x), y = quote(data$y), col = quote(scales::alpha("black", alpha))), moreargs))
    } else {
      do.call(plot, c(list(x = quote(data$x), y = quote(data$y), col = quote(scales::alpha(group_colors[data$group], alpha))), moreargs))
    }
    axis(3, labels = F, tck = 0.01)
    axis(4, labels = F, tck = 0.01)
    box()
    
    if(lm_show == TRUE & !is.null(lm_formula)){
      mapply(function(group_set, group_color){
        lm_tmp = lm(lm_formula, data = group_set)
        x_coords = seq(min(group_set$x), max(group_set$x), length.out = 100)
        y_coords = predict(lm_tmp, newdata = data.frame(x = x_coords))
        lines(x = x_coords, y = y_coords, col = group_color, lwd = 2.5)
      }, data_split, rgb(t(ceiling(col2rgb(group_colors)*0.8)), maxColorValue = 255))
    }
    
    # right density plot
    par(mar = c(4,0.25,0,1))
    plot(NULL, type = "n", ylim = moreargs$ylim, xlim = c(0, max(sapply(data_split, function(group_set) max(density(group_set$y, bw = bw)$y)))), main = NA, axes = F, xlab = "density")
    mapply(function(group_set, group_color){lines(x = density(group_set$y, bw = bw, adjust = adjust)$y, y = density(group_set$y, bw = bw)$x, col = group_color, lwd = 2)}, data_split, group_colors)
    axis(1)
  }, finally = {
    par(orig_par)
  })
}

# marginal_plot(x = Sepal.Width, y = Sepal.Length, group = Species, data = iris, lm_show = T)

## ggraster--ggplotclass--栅格数据可视化-可用于可视化环境变量和SDM投影####
ggraster <- function(data,occ,mcp = FALSE){
  ## data为env形式：最好划定范围；
  ## occ仅包含经纬度信息；
  if(mcp == TRUE){
    rang <- mcp(occ)
    data <- raster::crop(data,rang)
  }else{
    data =data
  }
  data2 <- as.data.frame(data,xy=T)
  names(data2)[3] <- 'value'
  names(occ) <- c("longitude", "latitude")
  titll <- toupper(deparse(substitute(data)))
  
  ggplot() +
    geom_raster(data = data2, aes(x = x, y = y, fill = value)) +
    geom_point(data=occ, aes(x=longitude, y=latitude ),size= 5, col='green', cex=0.1) +
    coord_quickmap() +
    theme_bw() + 
    scale_fill_gradientn(colours=rev(terrain.colors(10)),
                         na.value = "#CCCCCC7f")+
    ggtitle(titll)
}
# ggraster(data,occ,mcp= TRUE)

## 第二种：对栅格数据分组后可视化：####
ggplotclass <- function(data){
  class <- summary(round(data,3)) %>% data.frame()
  ## 分别获取最小、第一四分卫，中值，第三四分卫，最大值；
  vacl <- class[,1][1:5] 
  ## 构造class:
  vacl2 <- c(vacl[1:2],1,vacl[2:3],2,vacl[3:4],3,vacl[4:5],4)
  
  ff <- reclassify(data2,vacl2)
  
  library(rasterVis)
  gplot(ff) +
    geom_raster(aes(fill = factor(value))) +
    coord_equal()+
    guides(color=guide_legend(title = deparse(substitute(data))))
}
# ggplotclass(data2)
# 更简单的方法：
# rasterVis::levelplot(rprob, zscaleLog = NULL, contour = TRUE, margin = FALSE, at = miat, 
#           colorkey = myColorkey)







###############################################################
############ 第四部分：生态建模快捷命令 #######################
###############################################################



## filtenv--基于变量相关性、主成分分析和刀切法来筛选建模变量 ####
filtenv <- function(route,occ,route_out){ 
  # route 为环境数据集合的分布位置：
  # occ 仅包含经纬度的地理分布信息 ；
  # route_out:输出路径，到对应文件夹即可：
  
  tt <- route
  envsall <-  stack(tt)
  
  allva <- raster::extract(envsall,x) %>% na.omit() %>% data.frame()
  library(ade4) ## dudi.pca()
  library(factoextra) ## 可视化
  library(FactoMineR)
  pcaall <- PCA(allva )
  route <- paste0(route_out,"/",deparse(substitute(x)),"/pca.pdf")
  pdf(route)
  print(fviz_pca_var(pcaall, col.var = "black"))
  ## 主轴贡献排序：
  print(fviz_eig(pcaall, addlabels = TRUE, ylim = c(0, 100)))
  ## 查看对主轴1和2贡献最高的变量：
  print(fviz_contrib(pcaall, choice = "var", axes = 1))
  print(fviz_contrib(pcaall, choice = "var", axes = 2))
  dev.off()
  all_cor <- cor(allva,method = "pearson")
  library(pheatmap)
  route2 <- paste0(route_out,"/",deparse(substitute(x)),"/heatmap.png")
  png(route2,width =3400, height = 2000,res =300)
  pheatmap(all_cor,cluster_cols = F,cluster_rows = F,display_numbers = T,
           main = "Correlation Analysis of Environmental Factors of Origin_30s")
  dev.off()
  library(SDMtune)
  ## 构建研究范围：
  occs <- occ
  mcp <- function (xy) {
    library(rgdal)
    xy <- as.data.frame(coordinates(xy))
    coords.t <- chull(xy[, 1], xy[, 2])
    xy.bord <- xy[coords.t, ]
    xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
  }
  
  
  bgExt3 <- mcp(occs) %>% rgeos::gBuffer(.,width=1) %>% mask(crop(envsall[[1]],.),.)
  bg_xy2 <- dismo::randomPoints(bgExt3,3000,p=occs)
  
  bg_all <- as.data.frame(bg_xy2 ,colnames(c("long","lat")))
  
  data_sa <- prepareSWD(species = "hs ",
                        p =occs, a = bg_all,
                        env = envsall)
  
  sa_model <- train(method = "Maxent", data = data_sa, fc = "lqp", reg = 1, iter = 500)
  library(zeallot)  # For unpacking assignment
  c(train, test) %<-% trainValTest(data_sa, test = 0.4, only_presence = TRUE, seed = 25)
  
  jk <- doJk(sa_model, metric = "auc", test = test)
  jk_sa <- jk[,c(1,3)] %>% .[order(.$Train_AUC_withonly, decreasing= T),]
  route3 <- paste0(route_out,"/",deparse(substitute(x)),"/jk.txt")
  write.table(jk_sa,route3)
  
}
# route = "c:/Users/admin/Desktop/bio19"
# occ <- read.csv("hh.csv") ## 仅包含经纬度信息：
# route_out = "c:/Users/admin/Desktop/env_filter"
# filt(route,filt(route, occ, route_out)occ,route_out)


## tune2-- 基于maxent/biomod2进行模型调参，用于maxent和biomod2建模 ####
tune2 <- function(route,occraw,route_out){ ## c为env,d为occ
  library(diomod2)
  ## route 为保存环境变量的路径 
  ## occraw 为仅包含分部信息的经纬度；
  ## route_out：调参后输出的路径
  library(tidyverse)
  env <- list.files(route,pattern = ".tif",full.names = T) %>%  stack(c)

  occf <- function(x,y){
    occcle <- raster::extract(y,x) %>% round(4) %>% 
      cbind(x) %>% na.omit() %>% data.frame()
    tt <- dim(occcle)[2]
    dd <- tt-1
    occ <- occcle[,dd:tt] %>% data.frame()
    return(occ)
  }
  occ <- occf(occraw,env)
  
  ### 根据研究区的分布面积确定随机点数量 
  mcp <- function (xy) {
    library(rgdal)
    xy <- as.data.frame(coordinates(xy))
    coords.t <- chull(xy[, 1], xy[, 2])
    xy.bord <- xy[coords.t, ]
    xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
  }
  rang1 <- mcp(occ) 
  rang2 <- rgeos::gBuffer(rang1,width = 0.1)
  area1 <- raster::area(rang2) * 0.000001
  
  ned <- function(yy){
    needle = 0
    x = yy
    if(x<100){
      needle = 1000
      
    }else if(x>100 && x<500){
      needle = 2000
      
    }else if(x>500 && x<1000){
      needle = 3000
      
    }else if(x>1000 && x<2000){
      needle = 4000
    }else{needle = 5000
    }
    return(needle)
  }
  ned_va = ned(area1)
  
  
  
  occ_pr <- SpatialPointsDataFrame(coords=occ, 
                                   data=occ, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  
  occ_data <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(occ_pr)),
    resp.xy = occ,
    expl.var = env,
    resp.name = "occall",
    PA.nb.rep = 1,
    PA.nb.absences = ned_va,
    PA.strategy = 'random',
    na.rm = TRUE
  )
  ## 模型调参：
  tt <- BIOMOD_tuning(
    occ_data,
    models = c("GBM", "ANN" ), # ,   "MARS", "RF"
    models.options = BIOMOD_ModelingOptions(),
    method.ANN = "avNNet",
    method.RF = "rf",
    method.MARS = "earth",
    trControl = NULL,
    metric = "ROC",
    ctrl.RF = NULL,
    ctrl.ANN = NULL,
    ctrl.MARS = NULL,
    ctrl.GBM = NULL,
    tuneLength = 30,
    decay.tune.ANN = c(0.001, 0.01, 0.05, 0.1),
    size.tune.ANN = c(2, 4, 6, 8),
    maxit.ANN = 500,
    MaxNWts.ANN = 10 * (ncol(occ_data@data.env.var) + 1) + 10 + 1,
    interaction.GLM = c(0, 1),
    cvmethod.ME = "randomkfold",
    overlap.ME = FALSE,
    kfolds.ME = 4,
    n.bg.ME = 5000,
    env.ME = NULL,
    metric.ME = "ROC",
    clamp.ME = TRUE,
    parallel.ME = TRUE,
    numCores.ME = 10,
    Yweights = NULL 
  )
  
  ## 生成随机分布点：
  ## 统计面积返回随机点的数量：
  ## 基于分布构建点构建mcp：
  
  
  
  bg_xy2 <- sp::spsample(rang2, n=ned_va, type="random") 
  bg_xy <-  as.data.frame(bg_xy2@coords)
  
  library(SDMtune)
  datasdm <- prepareSWD(species = "species", 
                        p =occ , a = bg.xy, 
                        env = env)
  
  ## 利用sdmTUNE构建最优参数选择：
  folds <- get.checkerboard1(occ = datasdm@coords[datasdm@pa == 1, ],
                             env = env,
                             bg.coords = datasdm@coords[datasdm@pa == 0, ],
                             aggregation.factor = 4)
  
  model <- train(method = "Maxnet", data = datasdm,folds = folds)
  
  h <- list(reg = seq(0.5, 4, 0.5), fc = c("lq", "lh","qph","lqp", "lqph")) # , 
  om <- optimizeModel(model, hypers = h, metric = "tss", seed = 4,gen=20,
                      pop = 10)
  max <- head(om@results)[1,]  # Best combinations
  
  aa <- list(tt$tune.GBM$bestTune,tt$tune.RF$bestTune,tt$tune.ANN$bestTune,tt$tune.MARS$bestTune,
             max)
  names(aa) <- c("GBM","RF","ann","mars","max")
  
  rout <- paste0(route_out,"/",deparse(substitute(c)),".csv")
  write.csv(aa,rout)
  return(aa)
  
}
## biomod3--运行biomod2 并返回建模投影及评估结果 ####
## 需要根据tune2返回 的结果填写下一步建模的参数：
# library(biomod2)
# opt2 <- 
#   BIOMOD_ModelingOptions(
#     GBM = list(n.trees = 350,interaction.depth = 4,
#                shrinkage=0.1 , n.minobsinnode=10  ),
#     RF = list(ntree =1000,mtry=3 ),
#     ANN = list(size = 8,deca=0.001 , bag= FALSE ),
#     MAXENT.Phillips = list(linear= TRUE,quadratic=TRUE,
#                            product=TRUE,hinge= TRUE,
#                            memory_allocated = 4096,
#                            betamultiplier = 2.5)
#   )
# xh_opt2 = opt2

biomod3 <- function(route,occraw,route_out,rep,opt){  ## 其中a表示环境数据集，b表示物种分布数据,c为opt;
  ## route 为保存环境变量的路径 
  ## occraw 为仅包含分部信息的经纬度；
  ## route_out：调参后输出的路径
  ## rep:模型重复运行的次数
  ## opt:模型运行所需要的参数
  opt2 = opt
  
  
  env <- stack(route)
  occf <- function(x,y){
    occcle <- raster::extract(y,x) %>% round(4) %>% 
      cbind(x) %>% na.omit() %>% data.frame()
    tt <- dim(occcle)[2]
    dd <- tt-1
    occ <- occcle[,dd:tt] %>% data.frame()
    return(occ)
  }
  occ <- occf(occraw,env)
  
  
  ### 根据研究区的分布面积确定随机点数量 
  mcp <- function (xy) {
    library(rgdal)
    xy <- as.data.frame(coordinates(xy))
    coords.t <- chull(xy[, 1], xy[, 2])
    xy.bord <- xy[coords.t, ]
    xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
  }
  rang1 <- mcp(occ) 
  rang2 <- rgeos::gBuffer(rang1,width = 0.1)
  area1 <- raster::area(rang2) * 0.000001
  
  ned <- function(yy){
    needle = 0
    x = yy
    if(x<100){
      needle = 1000
      
    }else if(x>100 && x<500){
      needle = 2000
      
    }else if(x>500 && x<1000){
      needle = 3000
      
    }else if(x>1000 && x<2000){
      needle = 4000
    }else{needle = 5000
    }
    return(needle)
  }
  ned_va = ned(area1)
  
  ### 根据env或者occ来构建随机姓名：
  nnn <- paste0(deparse(substitute(occraw)),c(1:3))
  

  occ_pr <- SpatialPointsDataFrame(coords=occ, 
                                   data=occ,
                                   proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  
  occ_data <- BIOMOD_FormatingData(
    resp.var = rep(1, nrow(occ_pr)),
    resp.xy = occ,
    expl.var = env,
    resp.name = nnn[1],
    PA.nb.rep = 1, ## 5
    PA.nb.absences = ned_va,
    PA.strategy = 'random',
    na.rm = TRUE
  )

  
  s_models <- 
    BIOMOD_Modeling(
      data = occ_data,
      models = c( "GBM", "RF","ANN","MAXENT.Phillips"),
      models.options = opt2 ,
      NbRunEval =rep, #NbRunEval	integer, number of Evaluation run.
      DataSplit = 70, ##数据分割，80%
      VarImport = 2,
      models.eval.meth = c('KAPPA', 'TSS', 'ROC'),
      modeling.id = nnn[2],
      do.full.models = FALSE)
  
  route1 <- paste0(route_out,"/",deparse(substitute(occraw)),"/single_eva_run_models.txt")
  sink(route1)
  get_evaluations(s_models)
  sink()
  ensemble_models <- 
    BIOMOD_EnsembleModeling(
      modeling.output = s_models,
      em.by = 'all',
      eval.metric = 'TSS',
      eval.metric.quality.threshold = 0.6,
      models.eval.meth = c('TSS','ROC','KAPPA'),
      prob.mean = TRUE,
      prob.cv = TRUE, 
      committee.averaging = FALSE,
      prob.mean.weight = TRUE,
      VarImport = 2
    )
  route5 <-  paste0(route_out,"/",deparse(substitute(occraw)),"/ense_eva_run_models.txt")
  sink(route5)
  get_evaluations(ensemble_models)
  sink()
  
  (myBiomodModelOut_var_import <- get_variables_importance(s_models))
  ## make the mean of variable importance by algorithm
  varimp <- apply(myBiomodModelOut_var_import, c(1,2), mean)
  route2 <-  paste0(route_out,"/",deparse(substitute(occraw)),"/single_varimp_models.csv")
  write.csv(varimp,route2)
  route3 <-  paste0(route_out,"/",deparse(substitute(occraw)),"/ense_varimp_models.csv")
  varimpense <- get_variables_importance(ensemble_models)
  write.csv(varimpense,route3)
  
  route4 <-  paste0(route_out,"/",deparse(substitute(occraw)),"/singe_eva_models.pdf")
  pdf(route4)
  FILTER1 <- models_scores_graph(  ##输出多模型比对结果可视化：
    s_models, 
    by = "models", 
    metrics = c("ROC","TSS"), 
    xlim = c(0.5,1), 
    ylim = c(0.5,1)
  )
  
  FILTER2 <- models_scores_graph( ##输出run1和run2以及mean的结果；
    s_models, 
    by = "cv_run" , 
    metrics = c("ROC","TSS"), 
    xlim = c(0.5,1), 
    ylim = c(0.5,1)
  )
  
  ## 查看整合建模的结果：
  FILTER3 <- models_scores_graph( ##输出run1和run2以及mean的结果；
    ensemble_models, 
    by = "models" , 
    metrics = c("ROC","TSS"), 
    xlim = c(0.5,1), 
    ylim = c(0.5,1)
  )
  dev.off()
  
  
  ensemble_models_proj <- 
    BIOMOD_EnsembleForecasting( EM.output = ensemble_models,
                                new.env = env,
                                proj.name= nn[3],
                                do.stack =TRUE )
  route6 <-  paste0(route_out,"/",deparse(substitute(occraw)),"/ense_pro.pdf")
  pdf(route6)
  plot(ensemble_models_proj)
  dev.off()
  t1 <-  paste0(route_out,"/",deparse(substitute(occraw)),"/ense_mean.tif")
  t2 <-  paste0(route_out,"/",deparse(substitute(occraw)),"/ense_cv.tif")
  t3 <-  paste0(route_out,"/",deparse(substitute(occraw)),"/ense_wmean.tif")
  
  writeRaster(ensemble_models_proj@proj@val[[1]],t1)
  writeRaster(ensemble_models_proj@proj@val[[2]],t2)
  writeRaster(ensemble_models_proj@proj@val[[3]],t3)
  
  
}
# route = "c:/Users/admin/Desktop/bio19"
# occraw <- read.csv("hh.csv") ## 仅包含经纬度信息：
# route_out = "c:/Users/admin/Desktop/env_filter"
# rep= 5 ## 一般为5就要计算很长时间了；
# opt = xh_opt2 ## 模型运行的参数；

# biomod3(route,occraw,route_out,rep,opt)


## mx2-- 运行maxent，并返回建模投影及评估结果#####
mx2 <- function(occraw,routenow,routenowfut,route_out){  
  ## routenow 为当前环境条件下的气候信息位置；
  ## routefut 为未来或者过去气候条件下的气候信息位置；
  ## occraw 为仅包含分部信息的经纬度；
  ## route_out：模型输出的路径

  
  occ <- occraw
  
  ### 根据研究区的分布面积确定随机点数量 
  mcp <- function (xy) {
    library(rgdal)
    xy <- as.data.frame(coordinates(xy))
    coords.t <- chull(xy[, 1], xy[, 2])
    xy.bord <- xy[coords.t, ]
    xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
    return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
  }
  rang1 <- mcp(occ) 
  rang2 <- rgeos::gBuffer(rang1,width = 0.1)
  area1 <- raster::area(rang2) * 0.000001
  
  ned <- function(yy){
    needle = 0
    x = yy
    if(x<100){
      needle = 1000
      
    }else if(x>100 && x<500){
      needle = 2000
      
    }else if(x>500 && x<1000){
      needle = 3000
      
    }else if(x>1000 && x<2000){
      needle = 4000
    }else{needle = 5000
    }
    return(needle)
  }
  ned_va = ned(area1)

  rang <- rang2
  
  envnow <- list.files(routenow,pattern = ".tif",full.names = T) %>% stack()
  
  env1 <- envnow %>% crop(.,rang) %>% mask(.,rang)
  
  envfut <- paste0(routefut,pattern = ".tif",full.names = T) %>% stack()
  
  env2 <- envfut %>% crop(.,rang) %>% mask(.,rang)
  
  names(env1) <- names(env2) <-  names(envnow) 
  
  fold <- kfold(occ, k=5)
  test <- occ[fold == 1, ]
  train <- occ[fold != 1, ]
  
  back <- dismo::randomPoints(env1[[1]],ned_va, p =occ)
  args <-  c('linear=false', 'product=true', 'quadratic=true', 'hinge=false',
             'threshold=false','betamultiplier= 0.5','randomseed=true','writebackgroundpredictions=true',
             "applythresholdrule=Maximum training sensitivity plus specificity","threads=10",
             'pictures=false','outputformat=logistic','replicates=1','replicatetype=Bootstrap')
  xm <- maxent(env1, train,a =back ,args=  args)
  
 
  route1 <- paste0(route_out,"/",deparse(substitute(occraw)),"_响应曲线.tif")
  tiff(file = route1, res = 300, width = 3000, height =3000, compression = "lzw")
  response(xm)
  dev.off()
  
  
  t <- dismo::evaluate(p=test , a=back,
                       model = xm,x=env1)
  route2 <- paste0(route_out,"/",deparse(substitute(occraw)),"_模型评估.txt")
  sink(route2 )
  t
  sink()
  
  options(java.parameters = "-Xmx10g" )
  ## 裁剪后建模：
  now <- dismo::predict(xm,env1)
  fut <- dismo::predict(xm,env2)
  
  route3 <- paste0(route_out,"/",deparse(substitute(occraw)),"_NOW投影.tif")
  route3 <- paste0(route_out,"/",deparse(substitute(occraw)),"_FUT投影.tif")
  writeRaster(now,route3)
  writeRaster(fut,route3)
}
## ecospat_linux--运行ecospat包，脚本运行，适合服务器 多核心 ####
## occs为一个分布数据的集合体，使用rbind函数构建，包含物种名和经纬度信息；
## envs指定环境文件的路径
## 设定缓冲区的梯度范围：
## cores：指定核数；
## route_out:建模结果的输出路径：
ecospat_linux <- function(occs,route_env,buf,cores,route_out){
  ## occs为一个分布数据的集合体，使用rbind函数构建，包含物种名和经纬度信息；
  ## envs指定环境文件的路径
  ## 设定缓冲区的梯度范围：
  ## cores：指定核数；
  ## route_out:建模结果的输出路径：
  library(raster)
  library(dismo)
  library(tidyverse)
  library(SDMtune)
  library(rgdal)
  library(ecospat)
  ## 加载物种分布数据 
  occs <- ooccs
  
  ## 加载环境宾信息:
  biox <- list.files(route_env,pattern = ".tif",full.names = T) %>% stack()
  env1 <- raster("./F2_ENVS/aseutif/bio1.tif")
  env2 <- raster("./F2_ENVS/aseutif/bio10.tif")
  
  envs <- stack(env1,env2)
  ## 分布数据重新修正：
  n1 <- table(occs[1]) %>% t() %>% as.data.frame() 
  n2 <- dim(n1)[1]
  n3 <- n1[2]$Var2 %>% as.character()  
  
  
  ### 获取实际发生点的值：
  n4 = list()
  for(i in 1:n2){
    n4[[i]] <- raster::extract(envs,occs[which(occs[1] == n3[i]),2:3]) %>% na.omit(.) %>% data.frame(.)
  }
  names(n4) <-  n3 
  
  ## 获取不同缓冲区尺度下的背景值：
  ## 将数据框转为列表：获取缓冲区条件下随机分布点的经纬度信息 
  ## 根据分布区的面积返回随机点的数量：
  ned_ret <- function(xy){
    occ <- xy
    
    mcp <- function (xy) {
      library(rgdal)
      xy <- as.data.frame(coordinates(xy))
      coords.t <- chull(xy[, 1], xy[, 2])
      xy.bord <- xy[coords.t, ]
      xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
      return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
    }
    
    rang1 <- mcp(occ) 
    rang2 <- rgeos::gBuffer(rang1,width = 0.1)
    area1 <- raster::area(rang2) * 0.000001
    ned <- function(yy){
      needle = 0
      x = yy
      if(x<100){
        needle = 1000
        
      }else if(x>100 && x<500){
        needle = 2000
        
      }else if(x>500 && x<1000){
        needle = 3000
        
      }else if(x>1000 && x<2000){
        needle = 4000
      }else{needle = 5000
      }
      return(needle)
    }
    ned_va = ned(area1)
    return(ned_va)
    
  }
  
  names(occs)[1] <- "group"
  n5 <- occs  %>% group_by(group) %>% group_split()
  
  buffer2 <- function(x,buffer, ned){
    occs <- x
    occs.xy <- occs[c('longitude', 'latitude')]
    sp::coordinates(occs.xy) <- ~ longitude + latitude
    bgExt <- occs.xy
    bgExt3 <- rgeos::gBuffer(bgExt,width = buffer)
    bg_xy2 <- sp::spsample(bgExt3, n= ned, type="regular")
    
    return(as.data.frame(bg_xy2@coords))
  }
  
  n6 <- list()
  
  for(i in 1:n2){
    n7 <- n5[[i]][,2:3]
    
    ned <- ned_ret(n7)
    
    buffer  = buf
    
    n6[[i]] <- buffer2(n7,buffer,ned = ned)
  }
  names(n6) <- n3 
  
  ## 获取背景值对应的信息值：
  n8 = list()
  for(i in 1:n2){
    n9 <- n6[[1]]
    n8[[i]] <- raster::extract(envs,n9) %>% na.omit(.) %>% data.frame(.)
  }
  names(n8) <-  paste0("bg",n3) 
  
  ### 构建0_25度的下pca环境背景体系;
  
  all.back.env <- do.call(rbind,n8)
  
  clim.bkg <- all.back.env
  
  pca.env <-dudi.pca(clim.bkg, center = T, scale = T, scannf = F, nf = 2)
  
  scores.bkg<- pca.env$li
  
  zlist=list()
  for(i in 1:n2){
    n10 <- n8[[1]]
    n11 <- n4[[1]]
    assign(paste0("scores.bkg",n3[i]),suprow(pca.env,n10)$lisup)
    assign(paste0("scores.occ",n3[i]),suprow(pca.env,n11)$lisup)
    assign(paste0("z",n3[i]),ecospat.grid.clim.dyn(scores.bkg,
                                                   get(paste0("scores.bkg",n3[i])),
                                                   get(paste0("scores.occ",n3[i])),
                                                   500))
    zlist[[i]] <- get(paste0("z",n3[i]))
  }
  
  ### 025鼠李沙棘生态位保守性测试：
  n.groups <- n2
  g.names <- n3
  g.codenames <- n3
  rep <- 100
  D025 <- matrix(nrow = n.groups, ncol = n.groups)
  equ025 <- matrix(nrow = n.groups, ncol = n.groups)
  sim025 <- matrix(nrow = n.groups, ncol = n.groups)
  rownames(D025) <- colnames(D025) <- g.codenames
  rownames(equ025) <- colnames(equ025) <- g.codenames
  rownames(sim025) <- colnames(sim025) <- g.codenames
  
  cores <-  cores
  
  for (i in 2:n.groups) {
    
    for (j in 1:(i - 1)) {
      
      x1 <- zlist[[i]]
      x2 <- zlist[[j]]
      
      # Niche overlap
      D025[i, j] <- ecospat.niche.overlap (x1, x2, cor = TRUE)$D
      
      # niche equ
      equ025[i, j] <- ecospat.niche.similarity.test (x1, x2, rep,
                                                     alternative = "greater",ncores = cores )$p.D
      equ025[j, i] <- ecospat.niche.similarity.test (x1, x2, rep,
                                                     alternative = "lower",ncores = cores )$p.D
      # Niche similarity 
      sim025[i, j] <- ecospat.niche.similarity.test (x1, x2, rep,
                                                     alternative = "greater",ncores = cores )$p.D
      sim025[j, i] <- ecospat.niche.similarity.test (x1, x2, rep,
                                                     alternative = "lower",ncores = cores )$p.D
      
    }
  }
  da2 <- date()
  da3 <- gsub(" ","",da2)
  
  write.csv(D025,paste0(route_out,"/",da3,"_D025.csv"))
  write.csv(equ025,paste0(route_out,"/",equda3,"_D025.csv"))
  write.csv(sim025,paste0(route_out,"/",simda3,"_D025.csv"))
}












## 












###############################################################
#################### 生态位 可视化 ############################
###############################################################
## pca_plot_2p--计算PCA轴的差异性--可视化生态位轴的差异：####
pca_plot_2p <- function(x,y){ ## n1表示原产地，n2表示入侵地；
  library(ade4)
  n1 <- x
  n2 <- y
  sa_as_sdm <- rbind(n1,n2) %>% data.frame()
  sa_as_pca <- dudi.pca(sa_as_sdm,center= T ,scale=TRUE, scannf = F, nf = 3)
  sa_as.nat <- sa_as_pca$li[1:dim(n1)[1],]
  sa_as.inv <- sa_as_pca$li[dim(n1)[1]+1:dim(sa_as_sdm)[1],]
  scnat1 <- sa_as.nat$Axis1  %>% round(4) 
  scnat1 <- scnat1[!is.na(scnat1)]
  scinv1 <- sa_as.inv$Axis1  %>% round(4)
  scinv1 <- scinv1[!is.na(scinv1)]
  ## 进行一步数据缩放处理：
  scores1 <- c(scnat1,scinv1) 
  scores1 <- (scores1-min(scores1))/(max(scores1)-min(scores1))
  
  scnat2 <- sa_as.nat$Axis2  %>% round(4) 
  scnat2 <- scnat2[!is.na(scnat2)]
  scinv2 <- sa_as.inv$Axis2  %>% round(4)
  scinv2 <- scinv2[!is.na(scinv2)]
  scores2 <- c(scnat2,scinv2)
  scores2 <- (scores2-min(scores2))/(max(scores2)-min(scores2))
  
  scnat3 <- sa_as.nat$Axis3  %>% round(4) 
  scnat3 <- scnat3[!is.na(scnat3)]
  scinv3 <- sa_as.inv$Axis3  %>% round(4)
  scinv3 <- scinv3[!is.na(scinv3)]
  scores3 <- c(scnat3,scinv3)
  scores3 <- (scores3-min(scores3))/(max(scores3)-min(scores3))
  
  class1<-c(rep("nat",length(scnat1)),rep("inv",length(scinv1)))
  pc1<- data.frame(scores1,class1)
  class2<-c(rep("nat",length(scnat2)),rep("inv",length(scinv2)))
  pc2<- data.frame(scores2,class2)
  class3<-c(rep("nat",length(scnat3)),rep("inv",length(scinv3)))
  pc3<- data.frame(scores3,class3)
  
  pc <- list()
  pc[[1]] <- pc1
  pc[[2]] <- pc2
  pc[[3]] <- pc3
  return(pc)
}
# rbind(hgyawv,hgyaev)
## 后续可视化：
# time1 <- pca_plot_2p(hgyawv,hgyaev) ## p1显著；
# pc1 <- time1 [[1]]
# pc1$class1 <- factor(pc1$class1,labels = c("AS","SA"))
# pc2 <- time1 [[2]]
# pc2$class2 <- factor(pc2$class2,labels = c("AS","SA"))
# pc3 <- time1 [[3]]
# pc3$class3 <- factor(pc2$class2,labels = c("AS","SA"))
# library(PMCMR)
# p1 <- posthoc.kruskal.nemenyi.test(scores1~class1,data=pc1, dist="Tukey")
# p2 <- posthoc.kruskal.nemenyi.test(scores2~class2,data=pc2, dist="Tukey")
# p3 <- posthoc.kruskal.nemenyi.test(scores3~class3,data=pc3, dist="Tukey")
# ## 主要差异是
# 
# ## 轴投影的可视化：
# ggplot(pc2,aes(x=scores2 ,fill=class2))+ geom_density(alpha =0.4)+
#   guides()+xlab("PC2")


## dyn--可视化优化ecospat包的二维生态位重叠：####
dyn <- function (z1, z2, quant, title = "", name.axis1 = "Axis 1", 
                 name.axis2 = "Axis 2", interest = 1, colz1 = "#00FF0050", 
                 colz2 = "#FF000050", colinter = "#0000FF50", 
                 colZ1 = "green3", colZ2 = "red3") 
{
  if (is.null(z1$y)) {
    R <- length(z1$x)
    x <- z1$x
    xx <- sort(rep(1:length(x), 2))
    y1 <- z1$z.uncor/max(z1$z.uncor)
    Y1 <- z1$Z/max(z1$Z)
    if (quant > 0) {
      Y1.quant <- quantile(z1$Z[which(z1$Z > 0)], probs = seq(0, 
                                                              1, quant))[2]/max(z1$Z)
    }
    else {
      Y1.quant <- 0
    }
    Y1.quant <- Y1 - Y1.quant
    Y1.quant[Y1.quant < 0] <- 0
    yy1 <- sort(rep(1:length(y1), 2))[-c(1:2, length(y1) * 
                                           2)]
    YY1 <- sort(rep(1:length(Y1), 2))[-c(1:2, length(Y1) * 
                                           2)]
    y2 <- z2$z.uncor/max(z2$z.uncor)
    Y2 <- z2$Z/max(z2$Z)
    if (quant > 0) {
      Y2.quant <- quantile(z2$Z[which(z2$Z > 0)], probs = seq(0, 
                                                              1, quant))[2]/max(z2$Z)
    }
    else {
      Y2.quant = 0
    }
    Y2.quant <- Y2 - Y2.quant
    Y2.quant[Y2.quant < 0] <- 0
    yy2 <- sort(rep(1:length(y2), 2))[-c(1:2, length(y2) * 
                                           2)]
    YY2 <- sort(rep(1:length(Y2), 2))[-c(1:2, length(Y2) * 
                                           2)]
    plot(x, y1, type = "n", xlab = name.axis1, ylab = "density of occurrence")
    polygon(x[xx], c(0, y1[yy1], 0, 0), col = colz1, border = 0)
    polygon(x[xx], c(0, y2[yy2], 0, 0), col = colz2, border = 0)
    polygon(x[xx], c(0, apply(cbind(y2[yy2], y1[yy1]), 1, 
                              min, na.exclude = TRUE), 0, 0), col = colinter, border = 0)
    lines(x[xx], c(0, Y2.quant[YY2], 0, 0), col = colZ2, 
          lty = "dashed")
    lines(x[xx], c(0, Y1.quant[YY1], 0, 0), col = colZ1, 
          lty = "dashed")
    lines(x[xx], c(0, Y2[YY2], 0, 0), col = colZ2)
    lines(x[xx], c(0, Y1[YY1], 0, 0), col = colZ1)
    segments(x0 = 0, y0 = 0, x1 = max(x[xx]), y1 = 0, col = "white")
    segments(x0 = 0, y0 = 0, x1 = 0, y1 = 1, col = "white")
    seg.cat <- function(inter, cat, col.unf, col.exp, col.stab) {
      if (inter[3] == 0) {
        my.col = 0
      }
      if (inter[3] == 1) {
        my.col = col.unf
      }
      if (inter[3] == 2) {
        my.col = col.stab
      }
      if (inter[3] == -1) {
        my.col = col.exp
      }
      segments(x0 = inter[1], y0 = -0.01, y1 = -0.01, x1 = inter[2], 
               col = my.col, lwd = 4, lty = 2)
    }
    cat <- ecospat.niche.dyn.index(z1, z2, intersection = quant)$dyn
    inter <- cbind(z1$x[-length(z1$x)], z1$x[-1], cat[-1])
    apply(inter, 1, seg.cat, col.unf = "#00FF0050", 
          col.exp = "#FF000050", col.stab = "#0000FF50")
  }
  if (!is.null(z1$y)) {
    z <- t(as.matrix(z1$w + 2 * z2$w))[, nrow(as.matrix(z1$z.uncor)):1]
    z1$Z <- t(as.matrix(z1$Z))[, nrow(as.matrix(z1$Z)):1]
    z2$Z <- t(as.matrix(z2$Z))[, nrow(as.matrix(z2$Z)):1]
    if (interest == 1) {
      col1 = colorRampPalette(c("white", "black"),bias=1,alpha=0.8)(n=100)
      col2 = colorRampPalette(c("white", "black"),bias=1)(n=100)
      
      image(x = z2$x, y = z2$y, z = t(as.matrix(z2$z.uncor))[, 
                                                             nrow(as.matrix(z2$z.uncor)):1], col = col2, 
            zlim = c(1e-05, cellStats(z2$z.uncor, "max")))
      image(x = z1$x, y = z1$y, z = t(as.matrix(z1$z.uncor))[, 
                                                             nrow(as.matrix(z1$z.uncor)):1], col = col1, 
            zlim = c(1e-05,  cellStats(z1$z.uncor, "max")), 
            xlab = name.axis1, ylab = name.axis2, add = TRUE)
      
    }
    title(title)
    contour(x = z1$x, y = z1$y, z1$Z, add = TRUE, levels = quantile(z1$Z[z1$Z > 
                                                                           0], c(0, quant)), drawlabels = FALSE, lty = c(1, 
                                                                                                                         2), col = colZ1,lwd=3)
    contour(x = z2$x, y = z2$y, z2$Z, add = TRUE, levels = quantile(z2$Z[z2$Z > 
                                                                           0], c(0, quant)), drawlabels = FALSE, lty = c(1, 
                                                                                                                         2), col = colZ2,lwd=3)
    contour(x = z2$x, y = z2$y, z, add = TRUE,col="gray40",lwd=1.5,drawlabels = FALSE)
    
  }
}

## arow--可视化质心转移方向；
arow <- function (sp1, sp2, clim1, clim2, col = "##FFFFFF00") 
{
  if (ncol(as.matrix(sp1)) == 2) {
    arrows(median(sp1[, 1]), median(sp1[, 2]), median(sp2[, 
                                                          1]), median(sp2[, 2]), col = "gray27", lwd = 4, 
           length = 0.1)
    arrows(median(clim1[, 1]), median(clim1[, 2]), median(clim2[, 
                                                                1]), median(clim2[, 2]), lty = "11", col = col, 
           lwd = 2, length = 0.1)
  }
  else {
    arrows(median(sp1), 0.025, median(sp2), 0.025, col = "black", 
           lwd = 4, length = 0.1)
    arrows(median(clim1), -0.025, median(clim2), -0.025, 
           lty = "11", col = col, lwd = 2, length = 0.1)
  }
}

## 可视化单个物种的生态位：
ecospat.plot.niche2 <- function (z, title = "", name.axis1 = "Axis 1", name.axis2 = "Axis 2", 
                                 cor = FALSE) 
{
  if (is.null(z$y)) {
    R <- length(z$x)
    x <- z$x
    xx <- sort(rep(1:length(x), 2))
    if (cor == FALSE) 
      y1 <- z$z.uncor/max(z$z.uncor)
    if (cor == TRUE) 
      y1 <- z$z.cor/max(z$z.cor)
    Y1 <- z$Z/max(z$Z)
    yy1 <- sort(rep(1:length(y1), 2))[-c(1:2, length(y1) * 
                                           2)]
    YY1 <- sort(rep(1:length(Y1), 2))[-c(1:2, length(Y1) * 
                                           2)]
    plot(x, y1, type = "n", xlab = name.axis1, ylab = "density of occurrence")
    polygon(x[xx], c(0, y1[yy1], 0, 0), col = "grey")
    lines(x[xx], c(0, Y1[YY1], 0, 0))
  }
  if (!is.null(z$y)) {
    col1 = colorRampPalette(c("white", "black"),bias=1,alpha=0.8)(n=100)
    col2 = colorRampPalette(c("white", ""),bias=1,alpha=0.1)(n=100)
    if (cor == FALSE) 
      image(x = z$x, y = z$y, z = t(as.matrix(z$z.uncor))[, 
                                                          nrow(as.matrix(z$z.uncor)):1], col = col2, 
            zlim = c(1e-06, cellStats(z$z.uncor, "max")), 
            xlab = name.axis1, ylab = name.axis2)
    if (cor == TRUE) 
      image(x = z$x, y = z$y, z = t(as.matrix(z$z.uncor))[, 
                                                          nrow(as.matrix(z$z.uncor)):1], col = gray(100:0/100), 
            zlim = c(1e-06, cellStats(z$z.cor, "max")), 
            xlab = name.axis1, ylab = name.axis2)
    z$Z <- t(as.matrix(z$Z))[, nrow(as.matrix(z$Z)):1]
    contour(x = z$x, y = z$y, z$Z, add = TRUE, levels = quantile(z$Z[z$Z > 
                                                                       0], c(0, 0.5)), 
            drawlabels = FALSE, lty = c(1, 2),col = "grey40")
  } ## #33FF99b2:浅绿； ## 浅蓝：steelblue1 ## 浅红 pink； ## 浅黄 khaki
  title(title)
}


## plot.niche.all2--可视化优化ecospat包的二维生态位，提供多物种生态位可视化####
plot.niche.all2 <- function(z, n.groups, g.names,
                            contornar = TRUE, 
                            densidade = TRUE,
                            quantis = 10,
                            back = TRUE, title = "",
                            g.colors, n = 5,
                            cor1) {  
  
  # Color func
  cor1 <- function(cores.i, n) {
    al <- seq(0,1,(1/n))
    cores <- numeric(length(n))
    for(i in 1:n) {    
      corespar <- col2rgb(cores.i)/255
      cores[i] <- rgb(corespar[1, ], corespar[2, ],
                      corespar[3, ], alpha = al[i])
    }
    return(cores)
  }
  
  
  a <- list() 
  for(i in 1:n.groups) {
    a[[i]] <- colorRampPalette(c("transparent", cor1(g.colors[i], n)),
                               alpha = TRUE)  
  }
  
  xlim <- c(min(sapply(z, function(x){min(x$x)})),
            max(sapply(z, function(x){max(x$x)})))
  
  ylim <- c(min(sapply(z, function(x){min(x$y)})),
            max(sapply(z, function(x){max(x$y)})))
  
  image(z[[1]]$x, z[[1]]$y, as.matrix(z[[1]]$z.uncor), col = "white", 
        ylim = ylim, xlim = xlim,
        zlim = c(0.000001, max(as.matrix(z[[1]]$Z), na.rm = T)), 
        xlab = "PC1", ylab = "PC2", cex.lab = 1.5,
        cex.axis = 1.4)
  abline(h = 0, v = 0, lty = 2)
  box()
  
  if (back) {
    for(i in 1:n.groups) {
      contour(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$Z), add = TRUE,
              levels = quantile(z[[i]]$Z[z[[i]]$Z > 0], c(0, 1)),
              drawlabels = FALSE,lty = c(1, 2),
              col =  g.colors[i], lwd = 1)
    }
  }
  
  if (densidade) {
    for(i in 1:n.groups) {
      image(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$z.uncor),
            col = a[[i]](100), add = TRUE)
    }
  }
  
  
  if(contornar){
    for(i in 1:n.groups) {
      contour(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$z.uncor), add = TRUE,
              levels = quantile(z[[i]]$z.uncor[z[[i]]$z.uncor > 0],
                                seq(0, 1, (1/quantis)))[quantis],
              drawlabels = FALSE, lty = rev(c(rep(2,(quantis - 1)), 1)),
              col = rev(cor1(g.colors[i], quantis)),
              lwd = rev(c(rep(1, (quantis - 1)), 2)))
    }
  }
  
}
## plot.niche.single--单物种可视化 
plot.niche.single <- function(z, name.axis1 = "PC1", name.axis2 = "PC2",
                              cor = F, corte,  contornar = TRUE, 
                              densidade = TRUE, quantis = 10, 
                              back = TRUE, x = "red", title = "",
                              i) {  
  
  
  cor1 <- function(cores.i, n) {
    al <- seq(0,1,(1/n))
    cores <- numeric(length(n))
    for(i in 1:n) {    
      corespar <- col2rgb(cores.i)/255
      cores[i] <- rgb(corespar[1, ], corespar[2, ],
                      corespar[3, ], alpha = al[i])
    }
    return(cores)
  }
  
  
  a1 <- colorRampPalette(c("transparent",cor1(x, quantis)), alpha = TRUE)  
  
  xlim <- c(min(sapply(z, function(x){min(x$x)})),
            max(sapply(z, function(x){max(x$x)})))
  
  ylim <- c(min(sapply(z, function(x){min(x$y)})),
            max(sapply(z, function(x){max(x$y)})))
  
  graphics::image(z[[1]]$x, z[[1]]$y, as.matrix(z[[1]]$z.uncor), col = "white", 
                  ylim = ylim, xlim = xlim,
                  zlim = c(0.000001, max(as.matrix(z[[1]]$z.uncor), na.rm = T)), 
                  xlab = "PC1", ylab = "PC2", cex.lab = 1.5,
                  cex.axis = 1.4)
  
  abline(h = 0, v = 0, lty = 2)
  
  if (back) {
    contour(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$Z),
            add = TRUE, levels = quantile(z[[i]]$Z[z[[i]]$Z > 0],
                                          c(0, 0.5)), drawlabels = FALSE,
            lty = c(1, 2), col = x, lwd = 1)
  }
  
  if (densidade) {
    image(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$z.uncor), col = a1(100), add = TRUE)
  }
  
  
  if(contornar){
    contour(z[[i]]$x, z[[i]]$y, as.matrix(z[[i]]$z.uncor), 
            add = TRUE, levels = quantile(z[[i]]$z.uncor[z[[i]]$z.uncor > 0],
                                          seq(0, 1, (1 / quantis))),
            drawlabels = FALSE, lty = c(rep(2,(quantis - 1)), 1), 
            col = cor1(x, quantis), lwd = c(rep(1, (quantis - 1)), 2))
  }
  
  title(title)
  box()
}

















































































