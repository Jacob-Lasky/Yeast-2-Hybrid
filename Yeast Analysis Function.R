YH_Analysis <<- function(Directory, Bait_ID, Number_of_Replicates, Number_of_Plates_per_Replicate,Number_of_Wells){
                if(missing(Directory)){
                        print('Need directory')
                }
                if(missing(Bait_ID)){
                print('Please indicate the names of the transcription factors used')
                print('YH_Analysis(Directory,c(First Bait,Second Bait), Number_of_Replicates)')
                print('Quotes are required around each seperate Bait')
                stop()
                }
                if(missing(Number_of_Replicates)){
                        print('Please indicate the number of replicates for each transcription factor')
                        print('YH_Analysis(Directory,Bait_ID, c(x, y, z))')
                        print('Do not use quotes around each number')
                        stop()
                }
                if(missing(Number_of_Plates_per_Replicate)){
                        print('Please indicate the number of plates for each transcription factor')
                        print('YH_Analysis(Directory,Bait_ID, Number of Replcates, c(x, y, z))')
                        print('Do not use quotes around each number')
                        stop()
                }
                if(missing(Number_of_Wells)){
                        print('Please indicate the number of wells for each transcription factor')
                        print('YH_Analysis(Directory,Bait_ID, Number_of_Replicates, Number_of_Plates_per_Replicate, c(x, y, z))')
                        print('Do not use quotes around each number')
                        stop()
                }
                if(length(Number_of_Replicates == 1)){
                        Number_of_Replicates <- rep(Number_of_Replicates[1], length(Bait_ID))
                }
                if(length(Number_of_Plates_per_Replicate) == 1 ){
                        Number_of_Plates_per_Replicate <- rep(Number_of_Plates_per_Replicate[1], length(Bait_ID))
                }
                if(length(Number_of_Wells) == 1){
                        Number_of_Wells <- rep(Number_of_Wells[1], length(Bait_ID))
                }
        
        setwd(Directory)
        
        print('Reading TXT Files')
        Path <- Directory
        The_TXTs <- dir(Path,pattern = paste('\\.txt$',sep = ''),full.names = T)
        Short_Txts <<- dir(Path,pattern = paste('\\.txt$',sep = ''),full.names = F)
        x <- 1
        while(x <= (length(The_TXTs)/2)){
                Data <- read.delim(The_TXTs[x], stringsAsFactors = FALSE, skip = 7)
                Data <- Data[,-(length(Data))]
                assign(Short_Txts[x], Data)
                x <- x + 1
        }
        while(x <= length(The_TXTs)){
                Data <- read.delim(The_TXTs[x], stringsAsFactors = FALSE, skip = 6)
                Data <- Data[,-(length(Data))]
                assign(Short_Txts[x], Data)
                x <- x + 1
        }
#Import_Data function takes raw data and outputs easily to understand raw data tables.
        Import_Data <<- function(){
                print('Begin Importing')
                Path <- Directory
                Robot_Outputs <- Short_Txts
                Data_Names <- vector()
                test <- vector('numeric')
                x <- 1
                t <- 1
                Bait <- 1
                DF <- 1
        #Takes data from 6 - 384 well plates and creates a data frame with each replicate.
                for(All_Data in Robot_Outputs){
                        Short_Txts_ <- Short_Txts[((length(Short_Txts)/2)+1):length(Short_Txts)]
                        y <- 1
                        while(y <= Number_of_Plates_per_Replicate[Bait]){
                                Data_Vector <- Short_Txts
                                Data_Vector_ <- Short_Txts_
                                Short_Name <- Short_Txts
                                Short_Name_ <- Short_Txts_
                                Sheet <- get(Data_Vector[DF])
                                Sheet_ <- get(Data_Vector_[DF])
                                For_OD600 <- as.matrix(Sheet_)
                                For_OD420 <- as.matrix(Sheet)
                        #UX
                                Columns_Data <- as.data.frame(For_OD600[,seq(1,ncol(For_OD600),2)])
                                Rows_Data <- as.vector(t((Columns_Data[seq(1,nrow(Columns_Data),2),])))
                                Columns_Data <- as.data.frame(For_OD420[,seq(1,ncol(For_OD420),2)])
                                Rows_Data1 <- as.vector(t(Columns_Data[seq(1,nrow(Columns_Data),2),]))
                                OD600_and_OD420 <- data.frame(Rows_Data,Rows_Data1)
                                names(OD600_and_OD420) <- c('OD600', 'OD420')
                                Data1 <- OD600_and_OD420
                        #UX+1
                                Columns_Data <- as.data.frame(For_OD600[,seq(2,ncol(For_OD600),2)])
                                Rows_Data <- as.vector(t(Columns_Data[seq(1,nrow(Columns_Data),2),]))
                                Columns_Data <- as.data.frame(For_OD420[,seq(2,ncol(For_OD420),2)])
                                Rows_Data1 <- as.vector(t(Columns_Data[seq(1,nrow(Columns_Data),2),]))
                                OD600_and_OD420 <- data.frame(Rows_Data,Rows_Data1)
                                names(OD600_and_OD420) <- c('OD600', 'OD420')
                                Data2 <- OD600_and_OD420
                        #UX+2
                                Columns_Data <- as.data.frame(For_OD600[,seq(1,ncol(For_OD600),2)])
                                Rows_Data <- as.vector(t(Columns_Data[seq(2,nrow(Columns_Data),2),]))
                                Columns_Data <- as.data.frame(For_OD420[,seq(1,ncol(For_OD420),2)])
                                Rows_Data1 <- as.vector(t(Columns_Data[seq(2,nrow(Columns_Data),2),]))
                                OD600_and_OD420 <- data.frame(Rows_Data,Rows_Data1)
                                names(OD600_and_OD420) <- c('OD600', 'OD420')
                                Data3 <- OD600_and_OD420
                        #UX+3
                                Columns_Data <- as.data.frame(For_OD600[,seq(2,ncol(For_OD600),2)])
                                Rows_Data <- as.vector(t(Columns_Data[seq(2,nrow(Columns_Data),2),]))
                                Columns_Data <- as.data.frame(For_OD420[,seq(2,ncol(For_OD420),2)])
                                Rows_Data1 <- as.vector(t(Columns_Data[seq(2,nrow(Columns_Data),2),]))
                                OD600_and_OD420 <- data.frame(Rows_Data,Rows_Data1)
                                names(OD600_and_OD420) <- c('OD600', 'OD420')
                                Data4 <- OD600_and_OD420
                                DataX_Combo <- rbind.data.frame(Data1,Data2,Data3,Data4)
                                assign(paste(Bait_ID[Bait],' Rep ',t,' Plate ',y,sep = ""), DataX_Combo)
                                Data_Names <- c(Data_Names,paste(' Rep',t,' Plate ',y,sep = ''))
                                y <- y+1
                                DF <- DF + 1
                        }
                        ifelse(x == length(Data_Vector), x <- 1, x <- x+1)
                        if(t == Number_of_Replicates[t]){
                                if(Bait == length(Bait_ID)){
                                        break
                                        }
                                else{Bait <- Bait + 1}
                                t <- 0
                        }
                        t <- t+1
                        Data_Names <- Data_Names
                }
        #Names the data that has been imported into something that is easily understood
                for(All_Bait in Bait_ID){
                        l <- 1
                        Plates <- Number_of_Plates_per_Replicate[l]
                        New_Names <- vector()
                        while(l < (length(Short_Txts)/(2*Plates))+1){
                                Bait <- Bait_ID[l]
                                Replicates <- Number_of_Replicates[l]
                                R <- 1
                                while(R < (Replicates + 1)){
                                        Rep_Number <- 1
                                        Name <- paste(Bait,"(", Rep_Number,")",sep = '')
                                        New_Names <- c(New_Names,Name)    
                                        R <- R + 1
                                }
                                l <- l+1
                        }
                }
                File <- dir(path = Directory, pattern = paste('\\Identification',sep = ''), full.names = T)
                Identification <- read.csv(file = File[1], header = FALSE, na.strings=c("", "NA"))
                AGI_ <- as.data.frame(Identification[,1])
                Family_ <- as.data.frame(Identification[,2])
                AGI_and_Family <- cbind2(AGI_,Family_)
                AGI_and_Family <- AGI_and_Family[rowSums(is.na(AGI_and_Family)) != ncol(AGI_and_Family),]
                colnames(AGI_and_Family) <- c('AGI','Family')
                New_Names <<- New_Names
                Bait_Num <- 1
                Strain_Names <- vector()
                for (Bait in Bait_ID) {
                        Combined_Data <- matrix(ncol = 2)
                        colnames(Combined_Data) <- c('OD420', 'OD600')
                        Bait_Name <- Bait_ID[Bait_Num]
                        Rep_Num <- 1
                        while(Rep_Num < Number_of_Replicates[Bait_Num] + 1){
                                Plate_Num <- 1
                                while(Plate_Num < Number_of_Plates_per_Replicate[Bait_Num] + 1){
                                        Data_Name <- paste(Bait_ID[Bait_Num],'Rep',Rep_Num,'Plate',Plate_Num)
                                        Combined_Data <- rbind(Combined_Data, get(Data_Name))
                                        Plate_Num <- Plate_Num + 1
                                }
                                Rep_Num <- Rep_Num + 1
                        }
                        Combined_Data <- Combined_Data[-1,]
                        Well_Nums <- Number_of_Wells[Bait_Num]
                        Combined_Data <- Combined_Data[1:Well_Nums,]
                        Combined_Data <- Combined_Data[,c(2,1)]
                        as.numeric(Combined_Data[,1])
                        Combined_Data <- cbind(AGI_and_Family, Combined_Data)
                        assign(New_Names[Bait_Num], Combined_Data, envir = globalenv())
                        Bait_Num <- Bait_Num + 1
                }
        }
#Filter_Data function performs the calculations and outputs a data table that contains high confidence interactions.
        Filter_Data <<- function(Number_of_SD_for_Low_Confidence = 3, Number_of_SD_for_High_Confidence = 6, Number_of_SD_for_LUC_Over_OD600 = 3,
                                 Number_of_SD_for_Lum_Minus_BKG = 2){
                        print('Running Data Filter')
                        Average_NEG_OD600 <- vector()
                        Average_NEG_OD420 <- vector()
                        SD_NEG_OD600 <- vector()
                        SD_NEG_OD420 <- vector()
                        HC_Cutoff_OD600 <- vector()
                        LC_Cutoff_OD600 <- vector()
                        A_Data_Frame <- get(New_Names[1])
                        AGI_and_Family <- A_Data_Frame[,-(3:4)]
                        #Background calculation and outlier analysis of NEG wells
                        for(A_Strain in New_Names){
                                NEG_AGI <- get(A_Strain)
                                NEG_AGI <- NEG_AGI[which(NEG_AGI$AGI == 'NEG'),]
                                #OD600
                                Med_OD600 <- as.numeric(median(NEG_AGI$OD600))
                                NEG_OD600 <- as.numeric(unlist(NEG_AGI$OD600))
                                Middle_Subtracted_OD600 <- vector()
                                for(Value in NEG_OD600){
                                        Step1 <- Med_OD600 - Value
                                        Step2 <- abs(Step1)
                                        Middle_Subtracted_OD600 <- c(Middle_Subtracted_OD600,Step2)
                                }
                                Med_of_Med_Sub_OD600 <- median(Middle_Subtracted_OD600)
                                Non_Outliers_OD600 <- vector()
                                for(Value in NEG_OD600){
                                        if(abs(Med_OD600 - Value) <= (3.5*Med_of_Med_Sub_OD600)){
                                                Non_Outliers_OD600 <- c(Non_Outliers_OD600,Value)
                                        }
                                }
                                AVG_of_Non_Outliers_OD600 <- mean(Non_Outliers_OD600)
                                SD_of_Non_Outliers_OD600 <- sd(Non_Outliers_OD600)
                                Low_Confidence_Cutoff_OD600 <- (SD_of_Non_Outliers_OD600*Number_of_SD_for_Low_Confidence) + AVG_of_Non_Outliers_OD600
                                High_Confidence_Cutoff_OD600 <- (SD_of_Non_Outliers_OD600*Number_of_SD_for_High_Confidence) + AVG_of_Non_Outliers_OD600
                                #OD420
                                Med_OD420 <- as.numeric(median(NEG_AGI$OD420))
                                NEG_OD420 <- as.numeric(unlist(NEG_AGI$OD420))
                                Middle_Subtracted_OD420 <- vector()
                                for(Value in NEG_OD420){
                                        Step1 <- Med_OD420 - Value
                                        Step2 <- abs(Step1)
                                        Middle_Subtracted_OD420 <- c(Middle_Subtracted_OD420,Step2)
                                }
                                Med_of_Med_Sub_OD420 <- median(Middle_Subtracted_OD420)
                                Non_Outliers_OD420 <- vector()
                                for(Value in NEG_OD420){
                                        if(abs(Med_OD420 - Value) <= (3.5*Med_of_Med_Sub_OD420)){
                                                Non_Outliers_OD420 <- c(Non_Outliers_OD420,Value)
                                        }
                                }
                                AVG_of_Non_Outliers_OD420 <- mean(Non_Outliers_OD420)
                                SD_of_Non_Outliers_OD420 <- sd(Non_Outliers_OD420)
                                Average_NEG_OD600 <- c(Average_NEG_OD600,AVG_of_Non_Outliers_OD600)
                                Average_NEG_OD420 <- c(Average_NEG_OD420,AVG_of_Non_Outliers_OD420)
                                SD_NEG_OD600 <- c(SD_NEG_OD600,SD_of_Non_Outliers_OD600)
                                SD_NEG_OD420 <- c(SD_NEG_OD420,SD_of_Non_Outliers_OD420)
                                HC_Cutoff_OD600 <- c(HC_Cutoff_OD600,High_Confidence_Cutoff_OD600)
                                LC_Cutoff_OD600 <- c(LC_Cutoff_OD600,Low_Confidence_Cutoff_OD600)
                        }
                        AverageNEG_OD600 <- Average_NEG_OD600
                        AverageNEG_OD420 <- Average_NEG_OD420
                        SDNEG_OD600 <- SD_NEG_OD600
                        SDNEG_OD420 <- SD_NEG_OD420
                        High_C_Cutoff_OD600 <- HC_Cutoff_OD600
                        Low_C_Cutoff_OD600 <- LC_Cutoff_OD600
                        NEG_Data <- as.matrix(rbind(AverageNEG_OD420,AverageNEG_OD600,SDNEG_OD600,SDNEG_OD420,High_C_Cutoff_OD600,Low_C_Cutoff_OD600))
                        colnames(NEG_Data) <- New_Names
                        Cutoff_Values <<- NEG_Data
                        #Finding high confidence interactions
                        x <- 1
                        Cutoffs1 <- vector()
                        Cutoffs2 <- vector()
                        for(B_Strain in New_Names){
                                The_Data <- get(B_Strain)
                                OD420 <- as.numeric(unlist(The_Data[,4]))
                                LumBKG <- Cutoff_Values[1,x]
                                Lum_Minus_BKG <- OD420 - LumBKG
                                Lum_Minus_BKGDF <- as.data.frame(OD420 - LumBKG)
                                Quartiles <- quantile(Lum_Minus_BKG)
                                Third_Quarter <- Quartiles[4]
                                First_Quarter <- Quartiles[2]
                                Upper_Bound <- as.numeric(Third_Quarter + (1.5*(Third_Quarter - First_Quarter)))
                                Lower_Bound <- as.numeric(First_Quarter - (1.5*(Third_Quarter - First_Quarter)))
                                Outliers <- data.frame()
                                Test <- vector()
                                for(Each_Value in Lum_Minus_BKG){
                                        New_Value <- ifelse(Each_Value >= Upper_Bound, NA, Each_Value)
                                        Final_Value <- ifelse(New_Value < Lower_Bound, 'Outler-DN', New_Value)
                                        Outliers <- rbind(Outliers,Final_Value)
                                }
                                Outliers[is.na(Outliers)] <- 'Outlier-UP'
                                Outlier_OD420 <- Outliers
                                OD600 <- as.numeric(unlist(The_Data[,3]))
                                OD600DF <- as.data.frame(as.numeric(unlist(The_Data[,3])))
                                OD600BKG <- Cutoff_Values[2,x]
                                OD600_Minus_BKG <- as.data.frame(OD600 - OD600BKG)
                                OD600_Without_Outliers <- as.data.frame(ifelse(OD600DF[,1] < Cutoff_Values[6,x],NA,OD600_Minus_BKG[,1]))
                                Ratio_OD420andOD600 <- Lum_Minus_BKG / OD600_Minus_BKG
                                Analysis_Data <- data.frame(Lum_Minus_BKGDF,Outlier_OD420,OD600_Minus_BKG,OD600_Without_Outliers,Ratio_OD420andOD600)
                                No_Growth_Analysis <- as.data.frame(ifelse(is.na(Analysis_Data[,4]) == TRUE,NA,Analysis_Data[,5]))
                                Analysis_Data <- cbind(Analysis_Data,No_Growth_Analysis)
                                Quartiles <- quantile(unlist(Ratio_OD420andOD600))
                                Third_Quarter <- Quartiles[4]
                                First_Quarter <- Quartiles[2]
                                Upper_Bound <- as.numeric(Third_Quarter + (1.5*(Third_Quarter - First_Quarter)))
                                Lower_Bound <- as.numeric(First_Quarter - (1.5*(Third_Quarter - First_Quarter)))
                                Outliers <- data.frame()
                                for(Each_Value in Ratio_OD420andOD600){
                                        New_Value <- ifelse(Each_Value > Upper_Bound, NA, Each_Value)
                                        Final_Value <- ifelse(is.na(New_Value) == TRUE, 'Outlier-UP',ifelse(New_Value < Lower_Bound, 'Outler-DN', New_Value))
                                        Outliers <- rbind(Outliers,Final_Value)
                                }
                                Outlier_Ratios <- t(Outliers)
                                Analysis_Data <- cbind(Analysis_Data,Outlier_Ratios)
                                Numeric_Ratio_Only <- suppressWarnings(as.numeric(as.character(Outlier_Ratios)))
                                Ratio_Mean <- mean(Numeric_Ratio_Only, na.rm = TRUE)
                                Ratio_SD <- sd(Numeric_Ratio_Only, na.rm = TRUE)
                                Ratio_Cutoff <- Ratio_Mean + (Ratio_SD*Number_of_SD_for_LUC_Over_OD600)
                                Cutoffs1 <- c(Cutoffs1, Ratio_Cutoff)
                                Numeric_Lum_Minus_BKG_Only <- suppressWarnings(as.numeric(as.character(Analysis_Data[,2])))
                                LMB_Mean <- mean(Numeric_Lum_Minus_BKG_Only, na.rm = TRUE)
                                LMB_SD <- sd(Numeric_Lum_Minus_BKG_Only, na.rm = TRUE)
                                LMB_Cutoff <- LMB_Mean + (LMB_SD*Number_of_SD_for_Lum_Minus_BKG)
                                Cutoffs2 <- c(Cutoffs2, LMB_Cutoff)
                                #Slight rounding difference in means, standard deviations, and cutoffs than in XLSX
                                Ratio_Over_RatioCutoff <- unlist(Analysis_Data[,5]) / Ratio_Cutoff
                                FoldOver_Cutoff <- ifelse(is.na(Analysis_Data[,6]) == TRUE, NA, Ratio_Over_RatioCutoff)
                                Analysis_Data <- cbind(Analysis_Data, FoldOver_Cutoff)
                                x <- x + 1
                                #NA = No Growth
                                Likely_True_Positives <-       ifelse(is.na(Analysis_Data[,4]) == TRUE, NA, 
                                                                      ifelse(Analysis_Data[,8] <= .99, 'Negative',
                                                                             ifelse(Analysis_Data[,1] <= FoldOver_Cutoff, 'Potential False Positive',
                                                                                    'Likely True Positive')))
                                Analysis_Data <- cbind(AGI_and_Family, Analysis_Data, Likely_True_Positives)
                                Analysis_Data <- Analysis_Data[,-7]
                                colnames(Analysis_Data) <- c('AGI','Family','Lum-bkg','Lum-bkg (outliers)','OD600-bkg','Growth Analysis',
                                                             'Ratio Luc/OD600','Ratio (outliers)','Fold over Cutoff','Likely True Positives')
                                Analysis_Data[,8][is.na(Analysis_Data[,7])] <- NA
                                assign(B_Strain,Analysis_Data,envir=globalenv())
                                print(paste(B_Strain,'Filtered'))
                        }
                        Cutoffs1 <- t(data.frame(Cutoffs1))
                        rownames(Cutoffs1) <- 'Luc/OD600 (No Outliers)'
                        Cutoffs2 <- t(data.frame(Cutoffs2))
                        rownames(Cutoffs2) <- 'Lum-bkg (No Outliers)'
                        Cutoff_Values <<- rbind(Cutoff_Values, Cutoffs1, Cutoffs2)
        }
        Replicate_Analysis <<- function(){
                print('Running Replicate Analysis')
                Name_Pattern <- vector()
                The_Combined_Names <- vector()
                for(C_Strain in New_Names){
                        Pattern <- substr(C_Strain, 1, nchar(C_Strain)-2)
                        Pattern_Alt <- substr(C_Strain, 1, nchar(C_Strain)-3)
                        The_Combined_Names <- c(The_Combined_Names, paste('Combined-', Pattern_Alt, sep = ''))
                        Name_Pattern <- c(Name_Pattern,Pattern)
                }
                The_Combined_Names <<- The_Combined_Names
                All_Patterns <- Name_Pattern
                UniquePatterns <<- unique(All_Patterns)
                for(Each_Pattern in UniquePatterns){
                        Group <- New_Names[grep(Each_Pattern, New_Names, fixed = TRUE)]
                        assign(Each_Pattern, Group, envir = globalenv())
                }
                for(Each_Pattern in UniquePatterns){
                        Strain <- data.frame(y = 1:nrow(get(New_Names[1])))
                        for(Each_Group in get(Each_Pattern)){
                                The_DataMaster <- get(Each_Group)
                                The_DataFoldandPositives <- The_DataMaster[,9:10]
                                The_Data<- The_DataFoldandPositives
                                Strain <- cbind(Strain, The_Data)
                        }
                        The_DataAGIandFamily <- The_DataMaster[,1:2]
                        Strain <- Strain[,-1]
                        Strain <- cbind(The_DataAGIandFamily, Strain)
                        assign(Each_Pattern, Strain, envir = globalenv())
                }
                for(Each_Pattern in UniquePatterns){
                        Data <- get(Each_Pattern)
                        To_Maniuplate <- (Data[seq(from = 3, to = length(Data), by = 2)])
                #Shapiro-Wilks normality test
                        All_FoldOverCutoff <- unlist(To_Maniuplate)
                        All_FoldOverCutoff <- All_FoldOverCutoff[!is.na(All_FoldOverCutoff)]
                        if(length(All_FoldOverCutoff) > 5000){
                                All_FoldOverCutoff <- sample(All_FoldOverCutoff, 5000)
                        }
                        #print(shapiro.test(All_FoldOverCutoff))
                        High_Confidence <- Data
                #Plotting distribution
                        All_Y2H_Numbers <- density(All_FoldOverCutoff, na.rm = TRUE)
                        #plot(All_Y2H_Numbers)
                #Determining Mean, Stdev, and Coefficient of Variation
                        Mean <- rowMeans(To_Maniuplate, na.rm = TRUE)
                        Mean <- as.data.frame(Mean)
                        To_Maniuplate <- as.matrix(To_Maniuplate)
                        Standard_Deviation <- apply(To_Maniuplate,1,sd, na.rm = TRUE)
                        Standard_Deviation <- as.data.frame(Standard_Deviation)
                        Coefficient_of_Variation <- unlist(Mean) / unlist(Standard_Deviation)
                        Coefficient_of_Variation <- as.data.frame(Coefficient_of_Variation)
                        Combined <- cbind(Data, Mean, Standard_Deviation, Coefficient_of_Variation)
                        assign(Each_Pattern, Combined, envir = globalenv())
                }
                Number <- 1
                The_Combined_Names <- unique(The_Combined_Names)
                repeat{
                        if(is.na(UniquePatterns[Number]) == TRUE){break}
                        The_Data <- UniquePatterns[Number]
                        Data <- get(The_Data)
                        Name <- The_Combined_Names[Number]
                        assign(The_Combined_Names[Number], Data, envir = globalenv())
                        Number <- Number + 1
                }
                rm(list = UniquePatterns, envir = globalenv())
        }
        Export_Data <<- function(){
                x <- 1
                for(All_Files in New_Names){
                        Data <- get(New_Names[x])
                        write.csv(Data, paste(New_Names[x],'.csv', sep = ''), row.names = FALSE)
                        x <- x+1
                }
                print('Data Exported to CSV')
                print('NA Indicates No Growth')
        }

        
                
Import_Data()
        
Filter_Data()
        
#Replicate_Analysis()
        
Export_Data()
        
print('Analysis Complete')
}

YH_Analysis(Directory = "", 
            Bait_ID = c(""), 
            Number_of_Replicates = c(), 
            Number_of_Plates_per_Replicate = c(), 
            Number_of_Wells = c()
            )