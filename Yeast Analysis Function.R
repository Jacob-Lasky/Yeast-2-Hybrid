rm(list = ls())
YH_Analysis <<- function(Directory, Bait_ID, Number_of_Replicates, Number_of_Plates_per_Replicate,Number_of_Wells){
        x <- 1
        while(x == 1){
                if(missing(Directory)){
                        print('Need directory')
                        stop()
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
                x <- 0
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
                                For_Growth <- as.matrix(Sheet_)
                                For_Reporter <- as.matrix(Sheet)
                        #UX
                                Columns_Data <- as.data.frame(For_Growth[,seq(1,ncol(For_Growth),2)])
                                Rows_Data <- as.vector(t((Columns_Data[seq(1,nrow(Columns_Data),2),])))
                                Columns_Data <- as.data.frame(For_Reporter[,seq(1,ncol(For_Reporter),2)])
                                Rows_Data1 <- as.vector(t(Columns_Data[seq(1,nrow(Columns_Data),2),]))
                                Growth_and_Reporter <- data.frame(Rows_Data,Rows_Data1)
                                names(Growth_and_Reporter) <- c('Growth', 'Reporter')
                                Data1 <- Growth_and_Reporter
                        #UX+1
                                Columns_Data <- as.data.frame(For_Growth[,seq(2,ncol(For_Growth),2)])
                                Rows_Data <- as.vector(t(Columns_Data[seq(1,nrow(Columns_Data),2),]))
                                Columns_Data <- as.data.frame(For_Reporter[,seq(2,ncol(For_Reporter),2)])
                                Rows_Data1 <- as.vector(t(Columns_Data[seq(1,nrow(Columns_Data),2),]))
                                Growth_and_Reporter <- data.frame(Rows_Data,Rows_Data1)
                                names(Growth_and_Reporter) <- c('Growth', 'Reporter')
                                Data2 <- Growth_and_Reporter
                        #UX+2
                                Columns_Data <- as.data.frame(For_Growth[,seq(1,ncol(For_Growth),2)])
                                Rows_Data <- as.vector(t(Columns_Data[seq(2,nrow(Columns_Data),2),]))
                                Columns_Data <- as.data.frame(For_Reporter[,seq(1,ncol(For_Reporter),2)])
                                Rows_Data1 <- as.vector(t(Columns_Data[seq(2,nrow(Columns_Data),2),]))
                                Growth_and_Reporter <- data.frame(Rows_Data,Rows_Data1)
                                names(Growth_and_Reporter) <- c('Growth', 'Reporter')
                                Data3 <- Growth_and_Reporter
                        #UX+3
                                Columns_Data <- as.data.frame(For_Growth[,seq(2,ncol(For_Growth),2)])
                                Rows_Data <- as.vector(t(Columns_Data[seq(2,nrow(Columns_Data),2),]))
                                Columns_Data <- as.data.frame(For_Reporter[,seq(2,ncol(For_Reporter),2)])
                                Rows_Data1 <- as.vector(t(Columns_Data[seq(2,nrow(Columns_Data),2),]))
                                Growth_and_Reporter <- data.frame(Rows_Data,Rows_Data1)
                                names(Growth_and_Reporter) <- c('Growth', 'Reporter')
                                Data4 <- Growth_and_Reporter
                                DataX_Combo <- rbind.data.frame(Data1,Data2,Data3,Data4)
                                assign(paste(Bait_ID[Bait],' Rep ',t,' Plate ',y,sep = ""), DataX_Combo)
                                Data_Names <- c(Data_Names,paste(' Rep',t,' Plate ',y,sep = ''))
                                y <- y+1
                                DF <- DF + 1
                        }
                        ifelse(x == length(Data_Vector), x <- 1, x <- x+1)
                        if(t == Number_of_Replicates[grep(Bait_ID[Bait],Bait_ID)]){
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
                                Replicates <- Number_of_Replicates[grep(Bait,Bait_ID)]
                                R <- 1
                                while(R < (Replicates + 1)){
                                        Name <- paste(Bait,"(", R,")",sep = '')
                                        New_Names <- c(New_Names,Name)    
                                        R <- R + 1
                                }
                                if(l != length(Bait_ID)){
                                        l <- l+1
                                }
                                else{break}
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
                Rep_Num <- 1
                x <- 1
                while(Bait_Num <= length(Bait_ID)) {
                        Combined_Data <- matrix(ncol = 2)
                        colnames(Combined_Data) <- c('Reporter', 'Growth')
                        Bait_Name <- Bait_ID[Bait_Num]
                        Plate_Num <- 1
                        while(Plate_Num < Number_of_Plates_per_Replicate[Bait_Num] + 1){
                                Data_Name <- paste(Bait_ID[Bait_Num],'Rep',Rep_Num,'Plate',Plate_Num)
                                Combined_Data <- rbind(Combined_Data, get(Data_Name))
                                Plate_Num <- Plate_Num + 1
                        }
                        Combined_Data <- Combined_Data[-1,]
                        Well_Nums <- Number_of_Wells[Bait_Num]
                        Combined_Data <- Combined_Data[1:Well_Nums,]
                        Combined_Data <- Combined_Data[,c(2,1)]
                        as.numeric(Combined_Data[,1])
                        Combined_Data <- cbind(AGI_and_Family, Combined_Data)
                        assign(New_Names[x], Combined_Data, envir = globalenv())
                        assign(paste(New_Names[x],'-RAW',sep = ''), Combined_Data, envir = globalenv())
                        if(Rep_Num < Number_of_Replicates[Bait_Num]){
                                Rep_Num <- Rep_Num + 1
                        }
                        else{
                                Bait_Num <- Bait_Num + 1  
                                Rep_Num <- 1
                        }
                        x <- x+1
                }
        }
        #Filter_Data function performs the calculations and outputs a data table that contains high confidence interactions.
        Filter_Data <<- function(Number_of_SD_for_Low_Confidence = 3, Number_of_SD_for_High_Confidence = 6, Number_of_SD_for_Reporter_Over_Growth = 3,
                                 Number_of_SD_for_Reporter_Minus_BKG = 2){
                print('Running Data Filter')
                Average_NEG_Growth <- vector()
                Average_NEG_Reporter <- vector()
                SD_NEG_Growth <- vector()
                SD_NEG_Reporter <- vector()
                HC_Cutoff_Growth <- vector()
                LC_Cutoff_Growth <- vector()
                A_Data_Frame <- get(New_Names[1])
                AGI_and_Family <- A_Data_Frame[,-(3:4)]
                #Background calculation and outlier analysis of NEG wells
                for(A_Strain in New_Names){
                        NEG_AGI <- get(A_Strain)
                        NEG_AGI <- NEG_AGI[which(NEG_AGI$AGI == 'NEG'),]
                        #Growth
                        Med_Growth <- as.numeric(median(NEG_AGI$Growth))
                        NEG_Growth <- as.numeric(unlist(NEG_AGI$Growth))
                        Middle_Subtracted_Growth <- vector()
                        for(Value in NEG_Growth){
                                Step1 <- Med_Growth - Value
                                Step2 <- abs(Step1)
                                Middle_Subtracted_Growth <- c(Middle_Subtracted_Growth,Step2)
                        }
                        Med_of_Med_Sub_Growth <- median(Middle_Subtracted_Growth)
                        Non_Outliers_Growth <- vector()
                        for(Value in NEG_Growth){
                                if(abs(Med_Growth - Value) <= (3.5*Med_of_Med_Sub_Growth)){
                                        Non_Outliers_Growth <- c(Non_Outliers_Growth,Value)
                                }
                        }
                        AVG_of_Non_Outliers_Growth <- mean(Non_Outliers_Growth)
                        SD_of_Non_Outliers_Growth <- sd(Non_Outliers_Growth)
                        Low_Confidence_Cutoff_Growth <- (SD_of_Non_Outliers_Growth*Number_of_SD_for_Low_Confidence) + AVG_of_Non_Outliers_Growth
                        High_Confidence_Cutoff_Growth <- (SD_of_Non_Outliers_Growth*Number_of_SD_for_High_Confidence) + AVG_of_Non_Outliers_Growth
                        #Reporter
                        Med_Reporter <- as.numeric(median(NEG_AGI$Reporter))
                        NEG_Reporter <- as.numeric(unlist(NEG_AGI$Reporter))
                        Middle_Subtracted_Reporter <- vector()
                        for(Value in NEG_Reporter){
                                Step1 <- Med_Reporter - Value
                                Step2 <- abs(Step1)
                                Middle_Subtracted_Reporter <- c(Middle_Subtracted_Reporter,Step2)
                        }
                        Med_of_Med_Sub_Reporter <- median(Middle_Subtracted_Reporter)
                        Non_Outliers_Reporter <- vector()
                        for(Value in NEG_Reporter){
                                if(abs(Med_Reporter - Value) <= (3.5*Med_of_Med_Sub_Reporter)){
                                        Non_Outliers_Reporter <- c(Non_Outliers_Reporter,Value)
                                }
                        }
                        AVG_of_Non_Outliers_Reporter <- mean(Non_Outliers_Reporter)
                        SD_of_Non_Outliers_Reporter <- sd(Non_Outliers_Reporter)
                        Average_NEG_Growth <- c(Average_NEG_Growth,AVG_of_Non_Outliers_Growth)
                        Average_NEG_Reporter <- c(Average_NEG_Reporter,AVG_of_Non_Outliers_Reporter)
                        SD_NEG_Growth <- c(SD_NEG_Growth,SD_of_Non_Outliers_Growth)
                        SD_NEG_Reporter <- c(SD_NEG_Reporter,SD_of_Non_Outliers_Reporter)
                        HC_Cutoff_Growth <- c(HC_Cutoff_Growth,High_Confidence_Cutoff_Growth)
                        LC_Cutoff_Growth <- c(LC_Cutoff_Growth,Low_Confidence_Cutoff_Growth)
                }
                AverageNEG_Growth <- Average_NEG_Growth
                AverageNEG_Reporter <- Average_NEG_Reporter
                SDNEG_Growth <- SD_NEG_Growth
                SDNEG_Reporter <- SD_NEG_Reporter
                High_C_Cutoff_Growth <- HC_Cutoff_Growth
                Low_C_Cutoff_Growth <- LC_Cutoff_Growth
                NEG_Data <- as.matrix(rbind(AverageNEG_Reporter,AverageNEG_Growth,SDNEG_Growth,SDNEG_Reporter,High_C_Cutoff_Growth,Low_C_Cutoff_Growth))
                colnames(NEG_Data) <- New_Names
                Cutoff_Values <<- NEG_Data
                #Finding high confidence interactions
                x <- 1
                Cutoffs1 <- vector()
                Cutoffs2 <- vector()
                for(B_Strain in New_Names){
                        The_Data <- get(B_Strain)
                        Reporter <- as.numeric(unlist(The_Data[,4]))
                        ReporterBKG <- Cutoff_Values[1,x]
                        Reporter_Minus_BKG <- Reporter - ReporterBKG
                        Reporter_Minus_BKGDF <- as.data.frame(Reporter - ReporterBKG)
                        Quartiles <- quantile(Reporter_Minus_BKG)
                        Third_Quarter <- Quartiles[4]
                        First_Quarter <- Quartiles[2]
                        Upper_Bound <- as.numeric(Third_Quarter + (1.5*(Third_Quarter - First_Quarter)))
                        Lower_Bound <- as.numeric(First_Quarter - (1.5*(Third_Quarter - First_Quarter)))
                        Outliers <- data.frame()
                        Test <- vector()
                        for(Each_Value in Reporter_Minus_BKG){
                                New_Value <- ifelse(Each_Value >= Upper_Bound, NA, Each_Value)
                                Final_Value <- ifelse(New_Value < Lower_Bound, 'Outlier-DOWN', New_Value)
                                Outliers <- rbind(Outliers,Final_Value)
                        }
                        Outliers[is.na(Outliers)] <- 'Outlier-UP'
                        Outlier_Reporter <- Outliers
                        Growth <- as.numeric(unlist(The_Data[,3]))
                        GrowthDF <- as.data.frame(as.numeric(unlist(The_Data[,3])))
                        GrowthBKG <- Cutoff_Values[2,x]
                        Growth_Minus_BKG <- as.data.frame(Growth - GrowthBKG)
                        Growth_Without_Outliers <- as.data.frame(ifelse(GrowthDF[,1] < Cutoff_Values[6,x],NA,Growth_Minus_BKG[,1]))
                        Ratio_ReporterandGrowth <- Reporter_Minus_BKG / Growth_Minus_BKG
                        Analysis_Data <- data.frame(Reporter_Minus_BKGDF,Outlier_Reporter,Growth_Minus_BKG,Growth_Without_Outliers,Ratio_ReporterandGrowth)
                        No_Growth_Analysis <- as.data.frame(ifelse(is.na(Analysis_Data[,4]) == TRUE,NA,Analysis_Data[,5]))
                        Analysis_Data <- cbind(Analysis_Data,No_Growth_Analysis)
                        Quartiles <- quantile(unlist(Ratio_ReporterandGrowth))
                        Third_Quarter <- Quartiles[4]
                        First_Quarter <- Quartiles[2]
                        Upper_Bound <- as.numeric(Third_Quarter + (1.5*(Third_Quarter - First_Quarter)))
                        Lower_Bound <- as.numeric(First_Quarter - (1.5*(Third_Quarter - First_Quarter)))
                        Outliers <- data.frame()
                        for(Each_Value in Ratio_ReporterandGrowth){
                                New_Value <- ifelse(Each_Value > Upper_Bound, NA, Each_Value)
                                Final_Value <- ifelse(is.na(New_Value) == TRUE, 'Outlier-UP',ifelse(New_Value < Lower_Bound, 'Outlier-DOWN', New_Value))
                                Outliers <- rbind(Outliers,Final_Value)
                        }
                        Outlier_Ratios <- t(Outliers)
                        Analysis_Data <- cbind(Analysis_Data,Outlier_Ratios)
                        Numeric_Ratio_Only <- suppressWarnings(as.numeric(as.character(Outlier_Ratios)))
                        Ratio_Mean <- mean(Numeric_Ratio_Only, na.rm = TRUE)
                        Ratio_SD <- sd(Numeric_Ratio_Only, na.rm = TRUE)
                        Ratio_Cutoff <- Ratio_Mean + (Ratio_SD*Number_of_SD_for_Reporter_Over_Growth)
                        Cutoffs1 <- c(Cutoffs1, Ratio_Cutoff)
                        Numeric_Reporter_Minus_BKG_Only <- suppressWarnings(as.numeric(as.character(Analysis_Data[,2])))
                        LMB_Mean <- mean(Numeric_Reporter_Minus_BKG_Only, na.rm = TRUE)
                        LMB_SD <- sd(Numeric_Reporter_Minus_BKG_Only, na.rm = TRUE)
                        LMB_Cutoff <- LMB_Mean + (LMB_SD*Number_of_SD_for_Reporter_Minus_BKG)
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
                        colnames(Analysis_Data) <- c('AGI','Family','Reporter-bkg','Reporter-bkg (outliers)','Growth-bkg','Growth Analysis',
                                                     'Ratio Reporter/Growth','Ratio (outliers)','Fold over Cutoff','Likely True Positives')
                        Analysis_Data[,8][is.na(Analysis_Data[,7])] <- NA
                        assign(B_Strain,Analysis_Data,envir=globalenv())
                        print(paste(B_Strain,'Filtered'))
                }
                Cutoffs1 <- t(data.frame(Cutoffs1))
                rownames(Cutoffs1) <- 'Reporter/Growth (No Outliers)'
                Cutoffs2 <- t(data.frame(Cutoffs2))
                rownames(Cutoffs2) <- 'Reporter-bkg (No Outliers)'
                Cutoff_Values <<- rbind(Cutoff_Values, Cutoffs1, Cutoffs2)
        }
        Replicate_Analysis <<- function(){
        ###REALLY NEED TO WORK ON THIS
                if(all(Number_of_Replicates == 1)){
                        print('No Replicates Indicated')
                }
                else{
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
                        The_Combined_Names <<- unique(The_Combined_Names)
                        repeat{
                                if(is.na(UniquePatterns[Number]) == TRUE){break}
                                The_Data <- UniquePatterns[Number]
                                Data <- get(The_Data)
                                Name <- The_Combined_Names[Number]
                                assign(The_Combined_Names[Number], Data, envir = globalenv())
                                Number <- Number + 1
                        }
                        rm(list = UniquePatterns, envir = globalenv())
                        rm(UniquePatterns, envir = globalenv())
                        AVG_Growth <- vector()
                        AVG_Reporter <- vector()
                        AVG_Ratio <- vector()
                        for(Each_Name in New_Names){
                                Raw <- get(paste(Each_Name, "-RAW", sep = ''))
                                Data <- get(Each_Name)
                                RAW_AVG_Growth <- mean(Raw[,3])
                                RAW_AVG_Reporter <- mean(Raw[,4])
                                AVG_Ratio_ <- mean(Data[,7], na.rm = TRUE)
                                AVG_Growth <- c(AVG_Growth, RAW_AVG_Growth)
                                AVG_Reporter <- c(AVG_Reporter, RAW_AVG_Reporter)
                                AVG_Ratio <- c(AVG_Ratio, AVG_Ratio_)
                        }
                        AVGs <- rbind(AVG_Growth, AVG_Reporter, AVG_Ratio)
                        colnames(AVGs) <- New_Names
                        AVGs <<- AVGs
                        Total_Hits <- data.frame(nrow = Number_of_Wells[1])
                        for(Each_Name in New_Names){
                                Data <- get(Each_Name)
                                Hits <- unlist(Data[,10])
                                Total_Hits <- cbind(Total_Hits, Hits)
                        }
                        Total_Hits <- Total_Hits[,-1]
                        colnames(Total_Hits) <- New_Names
                        Total_Hits <- Total_Hits
                        x <- 1
                        Matching <- vector()
                        while(x <= length(Total_Hits)){
                                Positives <- match(Total_Hits[,x], 'Likely True Positive')
                                Matching <- cbind(Matching,Positives)
                                x <- x + 1
                        }
                        Matching <- Matching
                }
                
        }
        Export_Data <<- function(){
                Wells <- vector()
                U_Vals <- c(1:24)
                Letters_ <- LETTERS[1:8]
                Letter_Vals <- c(1:12)
                Vector_of_Well <- vector()
                x <- 1
                U_Num <- 1
                for(U in U_Vals){
                        for(LET in Letters_){
                                for(VALS in Letter_Vals){
                                        Vector_of_Well <- c(Vector_of_Well, paste("U", U, "-", LET, VALS, sep = ''))
                                }
                        }
                }
                Well_ID <- Vector_of_Well[1:Number_of_Wells[1]]
                x <- 1
                Control <- get(New_Names[1])
                Control <- Control[,9]
                y <- 1
                for(All_Files in New_Names){
                        Data <- get(New_Names[x])
                        Data_No_ID <- Data[,-(1:2)]
                        Data_Raw <- get(paste(New_Names[x],'-RAW',sep = ''))
                        Data <- cbind(Data_Raw, Data_No_ID)
                        Data <- cbind(Well_ID, Data, Control)
                        write.csv(Data, paste(New_Names[x],'.csv', sep = ''), row.names = FALSE)
                        x <- x + 1
                }
                if(all(Number_of_Replicates != 1)){
                        Data <- get(paste('Combined-', Bait_ID[y],sep = ''))
                        Data <- cbind(Well_ID, Data)
                        write.csv(Data, paste(The_Combined_Names[y],'.csv',sep = ''), row.names = FALSE)
                        y <- y + 1
                }
                print('Data Exported to CSV')
                print('NA Indicates No Growth')
        }
        
        
        
        Import_Data()
        
        Filter_Data()
        
        #Replicate_Analysis()
        
        Export_Data()
}
YH_Analysis(Directory = "D:/Documents/College/Research/SHINY TESTING/MCS 6", Bait_ID = c(1), Number_of_Replicates = c(1), Number_of_Plates_per_Replicate = c(6), Number_of_Wells = c(1999))
#Directory for CHE : "D:/Documents/College/Research/SHINY TESTING"
#Directory for MCS : "D:/Documents/College/Research/SHINY TESTING/MCS 6"
#Complete directory: 'D:/Documents/College/Research/Y2H/RAW/Uncompressed'

#Checklist
#       1. User directory is in quotes and is using forwardslash "/" not backslash "\"
#       2. Bait_ID values are indicated by 1:x, where x is the number of baits.
#       3. Number of replicates, plates, and wells are either all listed and seperated by commas
#               or, if this number is the same throughout, then listed a single time.
#       4a. Example: Bait_ID = c("Bait1"), Number_of_Replicates = c(2), Number_of_Plates_per_Replicate = c(6), Number_of_Wells = c(1999)
#               this will run "Bait1" as a duplicate and analyze it accordingly.
#       4b. Example: Bait_ID = c(1:6), Number_of_Replicates = c(1), Number_of_Plates_per_Replicate = c(6), Number_of_Wells = c(1999)
#       5. Make sure to include the Identification.csv file inside the folder that contains the TXT files that the robot gave
#       6. To easily run the program, check the box at the top left that says "Source on Save" and then press CTRL+S (CMD+S on Mac)
#       7. Make sure there is a .csv file in your directory that is named "Identification"

#Results will be given in .CSV files inside of the specified directory.


.libPaths('D:/Documents/College/Research/R-Programming/Rpackages')

String_Install <<- function(){
        source("http://bioconductor.org/biocLite.R")
        biocLite()
        source("http://bioconductor.org/biocLite.R")
        biocLite("STRINGdb")
}

if(length(grep("STRINGdb",library()))<1){String_Install()}


YH_Prioritization <<- function(BaitID, Bait_Number_to_Analyze, Directory){
        print('Begin Prioritization...please wait')
        New_Names <- New_Names[Bait_Number_to_Analyze]
        BaitID1 <- BaitID
        BaitID <- paste('3702.',BaitID,'.1',sep = '')
        library(STRINGdb)
        STRING_db <<- get_STRING_species(version = '10', species_name = 'Arabidopsis thaliana')
        STRING_db <<- STRINGdb$new( version="10", species=3702, score_threshold=0, input_directory="" )
        Sample <- get(New_Names[1])
        Sample <- as.character(Sample[,1])
        x <- 1
        STRING_ID <- vector()
        #while loop creates vector with STRING ID for each AGI
        while(x <= length(Sample)){
                NewID <- paste('3702.',Sample[x],'.1', sep = '')
                STRING_ID <- c(STRING_ID, NewID)
                x <- x + 1
        }
        #For loop adds the STRING IDs to each X(Y) set of data.
        STRING_ID <<- as.data.frame(STRING_ID)
        STRING_ID1 <- STRING_ID
        for(Each_Name in New_Names){
                Dat <- get(Each_Name)
                Dat <- cbind(STRING_ID, Dat)
                assign(Each_Name, Dat, envir = globalenv())
        }
        for(Each_Name in New_Names){
                POS <- get(Each_Name)
                POS <- as.character(POS[,length(POS)])
                Dat <- get(Each_Name)
                LIKELY_TRUE <- Dat[grep('Likely True Positive', POS),]
                assign(paste(Each_Name, 'Likely True Positives', sep = ''), LIKELY_TRUE, envir = globalenv())
        #Creates a csv file that shows the interactions between the "LIKELY TRUE POSITIVES" and bait
        #Creates a csv file that shows the interactions between all listed loci and bait
                if(file.exists(paste(Directory,'/', BaitID1,'_ALL_Interactions.csv',sep = '')) == FALSE){
                        Positives <- as.character(STRING_ID)
                        ALL_Interactions <- STRING_db$get_interactions(c(BaitID,Positives))
                        write.csv(ALL_Interactions, paste(Directory,'/', BaitID1,'_ALL_Interactions.csv',sep = ''))   
                }
                else{
                        ALL_Interactions <- read.csv(file = paste(Directory,'/', BaitID1,'_ALL_Interactions.csv',sep = ''), stringsAsFactors = FALSE)
                        #HC_Interactions <- ALL_Interactions
                }
                STRING_ALL_Interactions <- data.frame()
                BaitHit <- grep(BaitID, ALL_Interactions[,2])
                ALL_Interactions_FROM <- ALL_Interactions[BaitHit,]
                STRING_ALL_Interactions <- rbind(STRING_ALL_Interactions, ALL_Interactions_FROM)
                BaitHit <- grep(BaitID, ALL_Interactions[,3])
                ALL_Interactions_TO <- ALL_Interactions[BaitHit,]
                STRING_ALL_Interactions <- rbind(STRING_ALL_Interactions, ALL_Interactions_TO)
                STRING_ALL_Interactions <- STRING_ALL_Interactions[,-1]
                #HC_Interactions <- STRING_ALL_Interactions
                assign(paste('Bait_', Each_Name, '_ALL_Interactions', sep = ''), STRING_ALL_Interactions, envir = globalenv())
                #HC_Loci <- unique(unlist(STRING_HC_Interactions[,1:2]))
                HC_Loci <<- unique(unlist(STRING_ALL_Interactions[,1:2]))
                Combined_Scores <<- vector()
                STRING_COMBINED_SCORE <- STRING_ALL_Interactions
                for(Each_Loci in HC_Loci[-1]){
                        HC <- as.numeric(grep(Each_Loci, LIKELY_TRUE[,1]))
                        if(length(HC != 0)){
                                Combined_Score_Rows1 <- grep(Each_Loci, STRING_COMBINED_SCORE[,1])
                                Combined_Score_Rows2 <- grep(Each_Loci, STRING_COMBINED_SCORE[,2])
                                if(length(Combined_Score_Rows1) >= 1){
                                        Combined_Score_Rows <- Combined_Score_Rows1
                                }
                                if(length(Combined_Score_Rows2) >= 1){
                                        Combined_Score_Rows <- Combined_Score_Rows2
                                }
                                for(Each_Value in Combined_Score_Rows){
                                        Combined_Score <<- STRING_COMBINED_SCORE[Each_Value,16]
                                }
                                for(Each_Value in Combined_Score_Rows){
                                        STRING_COMBINED_SCORE <- STRING_COMBINED_SCORE[-(Each_Value),]
                                }
                                Combined_Scores <- c(Combined_Scores, Combined_Score)
                        }
                }
                Combined_Scores <<- Combined_Scores
                AVG_Combined_Score <- mean(Combined_Scores)
                AVG_Combined_Score <<- AVG_Combined_Score/1000
                rm(Combined_Scores)
                for(Each_Loci in HC_Loci){
                        HC <- as.numeric(grep(Each_Loci, LIKELY_TRUE[,1]))
                        if(length(HC) != 0){
                                HC_FoldoverCutoff <- LIKELY_TRUE[HC,10]
                                HC_Interactions <- get(paste('Bait_', Each_Name, '_ALL_Interactions', sep = ''))
                                #HC_Interactions <- get(paste(Each_Name, 'Likely True Positives', sep = ''))
                                HC_Interactions_LOC <- HC_Interactions[,1]
                                Combined_Score_LOC <- grep(Each_Loci, HC_Interactions_LOC)
                                if(length(Combined_Score_LOC) == 0){
                                        Combined_Score_LOC <- grep(Each_Loci, HC_Interactions[,2])
                                }
                                Combined_Score <- HC_Interactions[Combined_Score_LOC, 16]
                                Combined_Score_Ratio <- Combined_Score/1000
                                test1 <- Combined_Score_Ratio * 1000
                                Step1 <- test1 / AVG_Combined_Score
                                Adjusted_Yeast_Score <- HC_FoldoverCutoff^Step1
                                Combo <- cbind(Each_Loci, Combined_Score_Ratio, Adjusted_Yeast_Score)
                                Combo1 <- merge(Dat, Combo, by.x = 'STRING_ID', by.y = 'Each_Loci')
                                if(exists('HC_Combos') == TRUE){
                                        HC_Combos <- rbind(HC_Combos, Combo1)
                                }
                                else{
                                        HC_Combos <- Combo1
                                }
                        }
                }
                HC_Combos <- HC_Combos
                #Removes Bait-Bait self-interaction
                if(length(grep(BaitID, HC_Combos[,1])) != 0){
                        HC_Combos <- HC_Combos[-(grep(BaitID, x = HC_Combos[,1])),]
                }
                assign(paste('Bait_', Each_Name, '_HC_Combination', sep = ''), HC_Combos, envir = globalenv())
                write.csv(HC_Combos,paste(Directory, '/', 'High Confidence Adjusted Score.csv', sep = ''))
                rm(HC_Combos)                
                Combined_Scores <<- vector()
                rm(STRING_COMBINED_SCORE)
                rm(Combined_Score_Rows)
                rm(Combined_Score)
                STRING_COMBINED_SCORE <<- STRING_ALL_Interactions
                for(All_Loci in HC_Loci[-1]){
                                Combined_Score_Rows1 <- grep(All_Loci, STRING_COMBINED_SCORE[,1])
                                Combined_Score_Rows2 <- grep(All_Loci, STRING_COMBINED_SCORE[,2])
                                if(length(Combined_Score_Rows1) >= 1){
                                        Combined_Score_Rows <- Combined_Score_Rows1
                                }
                                if(length(Combined_Score_Rows2) >= 1){
                                        Combined_Score_Rows <- Combined_Score_Rows2
                                }
                                for(Each_Value in Combined_Score_Rows){
                                        Combined_Score <<- STRING_COMBINED_SCORE[Each_Value,16]
                                }
                                for(Each_Value in Combined_Score_Rows){
                                        STRING_COMBINED_SCORE <- STRING_COMBINED_SCORE[-(Each_Value),]
                                }
                                Combined_Scores <<- c(Combined_Scores, Combined_Score)
                }
                Combined_Scores <<- Combined_Scores
                #print(Combined_Scores)
                AVG_Combined_Score <- mean(Combined_Scores)
                AVG_Combined_Score <<- AVG_Combined_Score/1000
                All_Inters1 <- get(New_Names)
                All_Inters <- All_Inters1[,1]
                for(Each_Loci in HC_Loci){
                        All <- as.numeric(grep(Each_Loci, All_Inters))
                        if(length(All) != 0){
                                All_FoldoverCutoff <- as.numeric(All_Inters1[All,10])
                                All_Interactions <- get(paste('Bait_', Each_Name, '_ALL_Interactions', sep = ''))
                                #All_Interactions <- get(paste(Each_Name, 'Likely True Positives', sep = ''))
                                All_Interactions_LOC <- All_Interactions[,1]
                                Combined_Score_LOC <- grep(Each_Loci, All_Interactions_LOC)
                                if(length(Combined_Score_LOC) == 0){
                                        Combined_Score_LOC <- grep(Each_Loci, All_Interactions[,2])
                                }
                                Combined_Score <- All_Interactions[Combined_Score_LOC, 16]
                                Combined_Score_Ratio <- Combined_Score/1000
                                test1 <- Combined_Score_Ratio * 1000
                                Step1 <- test1 / AVG_Combined_Score
                                Adjusted_Yeast_Score <- All_FoldoverCutoff^Step1
                                Combo <- cbind(Each_Loci, Combined_Score_Ratio, Adjusted_Yeast_Score)
                                Combo1 <- merge(Dat, Combo, by.x = 'STRING_ID', by.y = 'Each_Loci')
                                if(exists('All_Combos') == TRUE){
                                        All_Combos <- rbind(All_Combos, Combo1)
                                }
                                else{
                                        All_Combos <- Combo1
                                }
                        }
                }
                All_Combos <- All_Combos
                #Removes Bait-Bait self-interaction
                if(length(grep(BaitID, All_Combos[,1])) != 0){
                        All_Combos <- All_Combos[-(grep(BaitID, x = All_Combos[,1])),]
                }
                assign(paste('Bait_', Each_Name, '_All_Combination', sep = ''), All_Combos, envir = globalenv())
                write.csv(All_Combos,paste(Directory, '/', 'Adjusted Score.csv', sep = ''))
                #Combines all data to have STRING analysis and non-STRING
                dat1 <- get(Each_Name)[1:6]
                dat2 <- get(Each_Name)[7:11]
                i <- sapply(dat2, is.factor)
                dat2[i] <- lapply(dat2[i], as.character)
                dat2 <- dat2
                dat <- cbind(dat1, dat2)
                Merged <- merge(All_Combos, dat, by = colnames(dat), all.y = TRUE)
                Numeric_Merged <- Merged[,12]
                Numeric_Merged1 <- as.numeric(levels(Numeric_Merged))[Numeric_Merged]
                Numeric_Merged <- Merged[,13]
                Numeric_Merged2 <- as.numeric(levels(Numeric_Merged))[Numeric_Merged]
                Merged <- Merged[,-(12:13)]
                Merged <- cbind(Merged, Numeric_Merged1, Numeric_Merged2)
                colnames(Merged) <- c('STRING_id', colnames(Merged[,2:11]), 'Adjusted_Yeast_Score', 'Combined_Score_Ratio')
        #Database retrieval
                Merged_IDs <- as.vector(Merged[,1])[1]
                #Pub_Med_Interactions_with_Bait <- list()
                #Pub_Med <- list()
                #for(Each_ID in Merged_IDs){
                #        Pub_Med_Interactions_with_Bait_ <- STRING_db$get_pubmed_interaction(Each_ID, BaitID1)
                #        Pub_Med_Interactions_with_Bait <- list(Pub_Med_Interactions_with_Bait, Pub_Med_Interactions_with_Bait_)
                #        Pub_Med_ <- STRING_db$get_pubmed(Each_ID)
                #        Pub_Med <- list(Pub_Med, Pub_Med_)
                #}
                #Pub_Med_Interactions_with_Bait <<- Pub_Med_Interactions_with_Bait[-1]
                #Pub_Med <<- Pub_Med[-1]
        #With pubmed, I want to have it so they can CLICK on the ID to get a list of interactions        
                
                #GO_ID <<- STRING_db$get_annotations()[,1:2]
        #With GO, CLICK on ID to get a list of GO IDs.
                
                
                Protein_Description <- STRING_db$get_proteins()
                colnames(Protein_Description) <- c('STRING_id',colnames(Protein_Description[2:4]))
                Merged <- merge(Protein_Description, Merged, by.x = 'STRING_id', all.y = TRUE)
                colnames(Merged) <- c(colnames(Merged[1:14]),'Adjusted Yeast Score', 'Combined Score Ratio')
                Name <<- paste(Each_Name,'_Complete', sep = '')
                Merged[,c(7,9,10,11,13,16)] <- round(Merged[,c(7,9,10,11,13,16)],3)
                Merged1 <- Merged[,1:7]
                Merged2 <- Merged[,9:11]
                Merged3 <- Merged[,13:16]
                To_Be_Mergedx <- as.vector(Merged[,8])
                To_Be_Merged1 <- data.frame(To_Be_Mergedx,stringsAsFactors = FALSE)
                x <- 1
                for(Each_Value in To_Be_Mergedx){
                        if((Each_Value != 'Outlier-DOWN') & Each_Value != "Outlier-UP"){
                                Each_Value <- round(as.numeric(Each_Value),3)
                                To_Be_Merged1[x,1] <- Each_Value
                        }
                        else{To_Be_Merged1[x,1] <- Each_Value}
                        x <- x+1
                }
                To_Be_Merged1 <- To_Be_Merged1
                colnames(To_Be_Merged1) <- 'Ratio (outliers)'
                x <- 1
                To_Be_Merged2 <- as.vector(Merged[,12])
                for(Each_Value in To_Be_Merged2){
                        if((Each_Value != 'Outlier-DOWN') & (Each_Value != "Outlier-UP") & (is.na(Each_Value) == FALSE)){
                                Each_Value <- round(as.numeric(Each_Value),3)
                                To_Be_Merged2[x] <- Each_Value
                        }
                        else{To_Be_Merged2[x] <- Each_Value}
                        x <- x+1
                }
                To_Be_Merged2 <- as.data.frame(To_Be_Merged2)
                colnames(To_Be_Merged2) <- 'Reporter-bkg (outliers)'
                Merged <<- cbind(Merged1, To_Be_Merged1, Merged2, To_Be_Merged2, Merged3)
                To_Be_Mergedy <- as.vector(Merged[,12])
                x <- 1
                for(Each_Value in To_Be_Mergedy){
                        if((Each_Value != 'Outlier-DOWN') & (Each_Value != "Outlier-UP") & (is.na(Each_Value) == FALSE)){
                                Each_Value <- round(as.numeric(Each_Value),3)
                                To_Be_Mergedy[x] <- Each_Value
                        }
                        else{To_Be_Mergedy[x] <- Each_Value}
                        x <- x+1
                }
                Merged[,12] <- To_Be_Mergedy
                assign(paste(Each_Name,'_Complete', sep = ''),Merged,envir = globalenv())
                rm(All_Combos)
                IDs <<- unique(as.character(Dat[,1]))
                ALL_Interactions <- read.csv(file = paste(Directory,'/', BaitID1,'_ALL_Interactions.csv',sep = ''), stringsAsFactors = FALSE)
                ALL_Interactions1 <- unique(as.character(ALL_Interactions[,2]))
                ALL_Interactions2 <- unique(as.character(ALL_Interactions[,3]))
                ALL_Interactions <- unique(c(ALL_Interactions1, ALL_Interactions2))
                Unmapped <- vector()
                for(Each_ID in IDs){
                        if(length(grep(Each_ID, ALL_Interactions)) == 0){
                                Unmapped <- c(Unmapped, Each_ID)
                        }
                }
                assign(paste('Bait_', Each_Name, '_Unmapped', sep = ''), Unmapped, envir = globalenv())
        }
        Pkg_Install <<- function(){
                source("https://bioconductor.org/biocLite.R")
                biocLite("ath1121501.db")
        }
        
        if(length(grep("ath1121501.db",library()))<1){Pkg_Install()}
        
        library(ath1121501.db)
        
        print("ATH1121501 Library Installed")
        
        AGI <- Merged[,5]
        
        MicroArray1 <- as.list(ath1121501ACCNUM[mappedkeys(ath1121501ACCNUM)])
        MatchingIndicies <- vector()
        for(Each_ID in AGI){
                MatchingIndicies <- append(MatchingIndicies, grep(Each_ID, MicroArray1))
        }
        MatchingIndicies <<- MatchingIndicies
        MicroArray <- MicroArray1[MatchingIndicies]
        
        MA_Names <- names(MicroArray)
        print("Microarray Success")
        
        Expression_Data <- vector()
        #load("D:/Documents/College/Research/ExpressionData/CopyOfAllNorm_6289Allgenes.Rdata")
        load("D:/Documents/College/Research/ExpressionData/CopyOfAllUnnorm_6289Allgenes.Rdata")
        
        #Norm_Data <- data.frame(matrix(NA, nrow = length(NormDatExpMatch[,MA_Names[1]]), ncol = 1))
        Unorm_Data <- data.frame(matrix(NA, nrow = length(UnormDatExpMatch[,MA_Names[1]]), ncol = 1))
        
        print("Beginning Correlated Expression Analysis")
        for(Each_Number in MA_Names){
                #Expression_Data <- as.data.frame(NormDatExpMatch[,Each_Number])
                Expression_Data <- as.data.frame(UnormDatExpMatch[,Each_Number])
                colnames(Expression_Data) <- MicroArray[Each_Number]
                #Norm_Data <- cbind(Norm_Data, Expression_Data)
                Unorm_Data <- cbind(Unorm_Data, Expression_Data)
        }
        #Norm_Data <- Norm_Data[,-1]
        Unorm_Data <- Unorm_Data[,-1]
        
        #Correlated_Expression <- cor(Norm_Data)
        Correlated_Expression <<- cor(Unorm_Data)
        
        NotInMA <- vector()
        for(Each_ID in AGI){
                Care <- grep(Each_ID, MicroArray)
                if(length(Care) == 0){
                        NotInMA <<- append(NotInMA, Each_ID)
                }
        }
        NotInMA <- NotInMA[which(NotInMA != "NEG")]
        NotInMA <- NotInMA[which(NotInMA != "pEXP-AD")]
        Bait_ID <- BaitID1
        MA_AGIs <- colnames(Correlated_Expression)
        Bait_Correlated_Expression <- vector()
        To_Merge2 <- Merged[,16]
        Expression_Correlation <- vector(length = length(To_Merge2))
        Merged <- cbind(Merged,Expression_Correlation)
        for(Each_AGI in MA_AGIs){
                Cor_Exp <- Correlated_Expression[which(colnames(Correlated_Expression) == Each_AGI),which(rownames(Correlated_Expression) == Bait_ID)]
                i <- which(Merged[,5] == Each_AGI)
                Merged[i,17] <- Cor_Exp
                Merged <- Merged
                Bait_Correlated_Expression <- append(Bait_Correlated_Expression, Cor_Exp)
        }
        colnames(Merged) <- c(colnames(Merged[,1:16]),paste("Correlated Expression with ",Bait_ID, sep = ""))
        Merged[,17] <- round(as.numeric(Merged[,17]),3)
        x <- 1
        for(Each_Value in Merged[,17]){
                if(Each_Value == 0){
                        Merged[x,17] <- NA
                }
                else{Merged[x,17] <- Each_Value}
                x <- x+1
        }
        Vals <- Merged[,8]
        x <- 1
        for(Each_Value in Vals){
                if(Each_Value != 'Outlier-UP'){
                        Vals[x] <- round(as.numeric(Each_Value),3)
                }
                else(Vals[x] <- Each_Value)
                x <- x + 1
        }
        Merged[,8] <- Vals
        Merged <<- Merged
        print('Prioritization Complete')
}
#YH_Prioritization(BaitID = 'AT5G08330', Bait_Number_to_Analyze = 1, Directory = "D:/Documents/College/Research/SHINY TESTING")
YH_Prioritization(BaitID = 'AT3G24500', Bait_Number_to_Analyze = 1, Directory = "D:/Documents/College/Research/SHINY TESTING/MCS 6")


