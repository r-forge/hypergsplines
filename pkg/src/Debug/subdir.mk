################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../backsolve.o \
../calculateModel.o \
../dataStructure.o \
../distributions.o \
../exhaustive.o \
../exhaustive2.o \
../functionWraps.o \
../getLogMargLikEst.o \
../getRho.o \
../glmExhaustive.o \
../glmGetSamples.o \
../glmLogMargLik.o \
../glmStochSearch.o \
../hyp2f1.o \
../iwls.o \
../linApproxDens.o \
../linalgInterface.o \
../logMargLik.o \
../random.o \
../stochSearch.o \
../sum.o \
../zdensity.o 

CPP_SRCS += \
../aggregateModelsTable.cpp \
../backsolve.cpp \
../calculateModel.cpp \
../dataStructure.cpp \
../distributions.cpp \
../exhaustive.cpp \
../exhaustive2.cpp \
../functionWraps.cpp \
../getLogMargLikEst.cpp \
../getRho.cpp \
../glmExhaustive.cpp \
../glmGetSamples.cpp \
../glmLogMargLik.cpp \
../glmStochSearch.cpp \
../hyp2f1.cpp \
../iwls.cpp \
../linApproxDens.cpp \
../linalgInterface.cpp \
../logMargLik.cpp \
../random.cpp \
../stochSearch.cpp \
../sum.cpp \
../zdensity.cpp 

OBJS += \
./aggregateModelsTable.o \
./backsolve.o \
./calculateModel.o \
./dataStructure.o \
./distributions.o \
./exhaustive.o \
./exhaustive2.o \
./functionWraps.o \
./getLogMargLikEst.o \
./getRho.o \
./glmExhaustive.o \
./glmGetSamples.o \
./glmLogMargLik.o \
./glmStochSearch.o \
./hyp2f1.o \
./iwls.o \
./linApproxDens.o \
./linalgInterface.o \
./logMargLik.o \
./random.o \
./stochSearch.o \
./sum.o \
./zdensity.o 

CPP_DEPS += \
./aggregateModelsTable.d \
./backsolve.d \
./calculateModel.d \
./dataStructure.d \
./distributions.d \
./exhaustive.d \
./exhaustive2.d \
./functionWraps.d \
./getLogMargLikEst.d \
./getRho.d \
./glmExhaustive.d \
./glmGetSamples.d \
./glmLogMargLik.d \
./glmStochSearch.d \
./hyp2f1.d \
./iwls.d \
./linApproxDens.d \
./linalgInterface.d \
./logMargLik.d \
./random.d \
./stochSearch.d \
./sum.d \
./zdensity.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/home/daniel/R/forge/hypergsplines/src" -I/usr/local/lib/R/include -I/home/daniel/R/library/RcppArmadillo/include -I/home/daniel/R/library/Rcpp/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


