################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../DriverBase.cpp \
../KVDriver.cpp \
../KVOR.cpp \
../KVORmod.cpp \
../TOV.cpp \
../Walecka.cpp \
../aux.cpp \
../eos.cpp \
../setconst.cpp \
../solve.cpp \
../wrap.cpp 

OBJS += \
./DriverBase.o \
./KVDriver.o \
./KVOR.o \
./KVORmod.o \
./TOV.o \
./Walecka.o \
./aux.o \
./eos.o \
./setconst.o \
./solve.o \
./wrap.o 

CPP_DEPS += \
./DriverBase.d \
./KVDriver.d \
./KVOR.d \
./KVORmod.d \
./TOV.d \
./Walecka.d \
./aux.d \
./eos.d \
./setconst.d \
./solve.d \
./wrap.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/const/levmar-2.6 -O0 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


