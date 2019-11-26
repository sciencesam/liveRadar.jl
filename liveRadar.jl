__precompile__()
#module liveRadar

println("Compiling Libraries ...")
using LibSerialPort
using MAT
using Plots
using DataStructures
using FFTW
using Statistics
using Suppressor



function liveDoppler(port,baud)
    # constants
    c = 299792458;      # (m/s) speed of light
    fs = 40000;         # (Hz) sampling frequency
    cpi = 0.50;         # (s) coherent processing interval 
    fc = 2440e6;        # (Hz) Center frequency (chirp is OFF)
    maxSpeed = 30;      # (m/s) maximum speed to display
    overlapFactor = 8;  # number of overlapped pulse windows (1 for no overlap)
    ovsDop = 4;         # oversample factor for Doppler axis
    # derived constants
    frameShiftTime = cpi/(overlapFactor-1); # time between frame starts
    N = round(Int,cpi*fs);  # number of samples per cpi/frame

    #N=800

    nfft = ovsDop*N; # number of samples in FFT
    maxSpeed = 500
    bufferSize = 300

    # initilize CircularBuffer
    b=preload_buffer(bufferSize,maxSpeed)
    show_waterfall(b)
    println("sleep")
    sleep(1)
    println("wake")
    # open serial port and enable serial DataStructures

    ser = open(port,baud)
    sp_set_dtr(ser.ref,SP_DTR_ON)
    println("Listening to $port @ $(baud)baud...")
    bytecount = 2*N÷overlapFactor # 2 bytes per short, N/ov samples in the window
    x=zeros(N)
    counter=0
    try
        while true
            while bytesavailable(ser) >=bytecount
                #println("s:")
                q,y=sp_blocking_read(ser.ref,bytecount,2)
                t=bytecombine.(y[1:2:end],y[2:2:end])
                if mean(t) > 7000
                    t=bytecombine.(y[2:2:end],y[1:2:end])
                end

                x=[x;t[t.<7000]]

                doppler!(b,x[end-N:end],nfft,maxSpeed)


                if counter%2==1
                    @sync show_waterfall(b)
                end
                counter=counter+1
            end
        end
    catch e
        # exit on CTRL-C
        if e isa InterruptException
            println("InterruptException")
            close(ser)
            println("closing")

            save_mat_file(x,"test.mat")
        else
            close(ser)
            println("meh")
            println(e)
        end      
    end  
end
function offlineDoppler(matfile)
    # -------------------- Setup constants and parameters --------------------
    c = 299792458;      # (m/s) speed of light
    pw = 20e-3;         # (s) pulse length
    fs = 40000;         # (Hz) sample rate
    cpi = 0.50;         # (s) coherent processing interval 
    ovsRng = 10;        # oversampling factor applied in range
    fStart = 2400e6;    # (Hz) LFM start frequency
    fStop = 2480e6;     # (Hz) LFM stop frequency
    nPulseCancel = 2;   # Number of pulses to use for canceller
    maxRange = 100;     # (m) maximum range to display
    overlapFactor = 8;  # number of overlapped pulse windows (1 for no overlap)
    ovsDop = 4;         # oversample factor for Doppler axis
    # ------------------------------------------------------------------------
    Np = round(Int,pw*fs)-1     # samples per pulse (-1 due to start flag)
    nfft = Np*ovsRng
    bw = fStop - fStop          # transmit bandwidth

    # read bytes as they come with live serial
    samples = read_mat(matfile)
    N = round(Int,cpi*fs);  # number of samples per cpi/frame
    bytecount = 2*N÷overlapFactor

    x=zeros(N)
    # initilize CircularBuffer
    maxSpeed = 500
    bufferSize = 300
    b=preload_buffer(bufferSize,maxSpeed)

    show_waterfall(b)
    sleep(1)
    current_index = 1 # arrays start at 1 in Julia
    for i=1:(bytecount÷2):round_to(length(samples),bytecount÷2)
        t=samples[i:i+bytecount÷2-1] # simulated byte read from serial
        #println(t)
        x=[x;t[t.<7000]] # simulate adding bytes to heap 
        # find pulses in the shorts since the last look
        doppler!(b,x[end-N:end],ovsDop*N,maxSpeed)
        show_waterfall(b)
    end
end
function liveRanging(port,baud)
    # -------------------- Setup constants and parameters --------------------
    c = 299792458;      # (m/s) speed of light
    pw = 20e-3;         # (s) pulse length
    fs = 40000;         # (Hz) sample rate
    cpi = 0.50;         # (s) coherent processing interval 
    ovsRng = 10;        # oversampling factor applied in range
    fStart = 2400e6;    # (Hz) LFM start frequency
    fStop = 2480e6;     # (Hz) LFM stop frequency
    nPulseCancel = 2;   # Number of pulses to use for canceller
    maxRange = 100;     # (m) maximum range to display
    overlapFactor = 8;  # number of overlapped pulse windows (1 for no overlap)
    # ------------------------------------------------------------------------
    Np = round(Int,pw*fs)-1     # samples per pulse (-1 due to start flag)
    #Np = Np÷2
    nfft = Np*ovsRng
    bw = fStop - fStop          # transmit bandwidth

    # read bytes as they come with live serial
    #samples = read_mat(matfile)
    N = round(Int,cpi*fs);  # number of samples per cpi/frame
    bytecount = 2*N÷overlapFactor

    x=zeros(N)
    # initilize CircularBuffer
    startS = 3996
    maxSpeed = 533
    bufferSize = 600
    b=preload_buffer(bufferSize,maxSpeed)
    previous_frame = zeros(Complex{Float64},size(b[1]))

    show_waterfall(b,true)
    sleep(1)
    current_index = 1 # arrays start at 1 in Julia
    counter=1


    ser = open(port,baud)
    sp_set_dtr(ser.ref,SP_DTR_ON)
    println("Listening to $port @ $(baud)baud...")
    
    try
        while true
            while bytesavailable(ser) >=bytecount
                #println("s:")
                q,y=sp_blocking_read(ser.ref,bytecount,2)
                t=bytecombine.(y[1:2:end],y[2:2:end])
                if mean(t) > 7000
                    t=bytecombine.(y[2:2:end],y[1:2:end])
                end

                x=[x;t[t.<7000]]

                pulses = findPulses(x,current_index,Np)
                if length(pulses) > 0
                    current_index = maximum(pulses)+1
                    for pulse in pulses
                        ranging!(b,x[pulse:pulse+Np-1],nfft,maxSpeed,1, previous_frame)
                        #show_waterfall(b,true)
                        if counter%15==0
                            @sync show_waterfall(b,false)
                        end
                        counter=counter+1
                    end   
                end
            end
        end
    catch e
        # exit on CTRL-C
        if e isa InterruptException
            println("InterruptException")
            close(ser)
            println("closing")

            save_mat_file(x,"test.mat")
        else
            close(ser)
            println("meh")
            println(e)
        end      
    end  
end
function simliveRanging(matfile)
    # -------------------- Setup constants and parameters --------------------
    c = 299792458;      # (m/s) speed of light
    pw = 20e-3;         # (s) pulse length
    fs = 40000;         # (Hz) sample rate
    cpi = 0.50;         # (s) coherent processing interval 
    ovsRng = 10;        # oversampling factor applied in range
    fStart = 2400e6;    # (Hz) LFM start frequency
    fStop = 2480e6;     # (Hz) LFM stop frequency
    nPulseCancel = 2;   # Number of pulses to use for canceller
    maxRange = 100;     # (m) maximum range to display
    overlapFactor = 8;  # number of overlapped pulse windows (1 for no overlap)
    # ------------------------------------------------------------------------
    Np = round(Int,pw*fs)-1     # samples per pulse (-1 due to start flag)
    println(Np)
    #Np = Np÷2
    nfft = Np*ovsRng
    bw = fStop - fStop          # transmit bandwidth

    # read bytes as they come with live serial
    samples = read_mat(matfile)
    N = round(Int,cpi*fs);  # number of samples per cpi/frame
    bytecount = 2*N÷overlapFactor
    x=[]
    #x=zeros(N)
    # initilize CircularBuffer
    startS = 3996
    maxSpeed = 533
    bufferSize = 300
    b=preload_buffer(bufferSize,maxSpeed)
    previous_frame = zeros(Complex{Float64},size(b[1]))

    show_waterfall(b,true)
    sleep(1)
    current_index = 1 # arrays start at 1 in Julia
    for i=1:(bytecount÷2):round_to(length(samples),bytecount÷2)
        t=samples[i:i+bytecount÷2] # simulated byte read from serial
        x=[x;t]
        #x=[x;t[t.<7000]] # simulate adding bytes to heap 
        # find pulses in the shorts since the last look
        pulses = findPulses(x,current_index,Np)
        if length(pulses) > 0
            current_index = maximum(pulses)+1

            for pulse in pulses
                #current_index = copy(pulse)+Np+1
                #print(current_index,' ')
                ranging!(b,x[pulse:pulse+Np-1],nfft,maxSpeed,startS, previous_frame)
                show_waterfall(b,true)

            end
        end
        #println(current_index)
    end
end
function round_to(x,n)
    return floor(Int,n*floor(x/n))
end
function testloop()
    x=[]
    samples = collect(1:99)
    bytecount = 8
    current_index=1
    for i=1:(bytecount÷2):round_to(length(samples),bytecount÷2)
        t=samples[i:i+bytecount÷2-1] # simulated byte read from serial
        println(t)
        x=[x;t]
        #x=[x;t[t.<7000]] # simulate adding bytes to heap 
        # find pulses in the shorts since the last look
        pulses = findPulses(x,current_index)
        for pulse in pulses
            if pulse+Np < length(x)
                println(x[pulse:pulse+Np-1])
                show_waterfall(b,true)
                current_index = pulse
            end
        end
    end
end
function offlineRanging(matfile)
    # -------------------- Setup constants and parameters --------------------
    c = 299792458;      # (m/s) speed of light
    pw = 20e-3;         # (s) pulse length
    fs = 40000;         # (Hz) sample rate
    cpi = 0.50;         # (s) coherent processing interval 
    ovsRng = 10;        # oversampling factor applied in range
    fStart = 2400e6;    # (Hz) LFM start frequency
    fStop = 2480e6;     # (Hz) LFM stop frequency
    nPulseCancel = 2;   # Number of pulses to use for canceller
    maxRange = 100;     # (m) maximum range to display
    overlapFactor = 8;  # number of overlapped pulse windows (1 for no overlap)
    # ------------------------------------------------------------------------
    Np = round(Int,pw*fs)-1     # samples per pulse (-1 due to start flag)
    nfft = Np*ovsRng
    bw = fStop - fStop          # transmit bandwidth

    # read bytes as they come with live serial
    samples = read_mat(matfile)
    N = round(Int,cpi*fs);  # number of samples per cpi/frame
    bytecount = 2*N÷overlapFactor

    x=zeros(N)
    # initilize CircularBuffer
    startS = 3996
    maxSpeed = 533
    bufferSize = 300
    b=preload_buffer(bufferSize,maxSpeed)
    previous_frame = zeros(Complex{Float64},size(b[1])) # copy b[1] for an array of zeros

    show_waterfall(b)
    sleep(1)
    current_index = 1 # arrays start at 1 in Julia
    x=samples[:]
    # find pulses in the shorts since the last look
    pulses = findPulses(x,current_index,Np)
    for pulse in pulses
        #if pulse+Np < length(x)
        ranging!(b,x[pulse:pulse+Np],nfft,maxSpeed,startS, previous_frame)
        show_waterfall(b,true)
        current_index = pulse
        #end
    end
end
function plotFrame(x)
    plot(x)
    gui()
end
function findPulses(x,i,n=0::Int)
    t=x[i:end]
    flag_pulsestart = 5000
    pulseStarts = t.>flag_pulsestart
    # findall returns indexes
    ind = findall(pulseStarts).+i # pulse start flags + current start index
    ind = ind[ind.+n.<length(x)]
    return ind
end
function bytecombine(firstByte,secondByte)
    return convert(Int32,convert(UInt16,firstByte)<<8+convert(UInt16,secondByte))
end
function hann_window(N)
    y = 0.5 .+ 0.5*cos.(2π*( (0:N-1)/(N-1) .- 0.5 ))
    return y
end
function mag(x)
    return 20*log10.(abs.(x))
end
function windowed_fft(x,nfft)
    # returns fft with hanning window applied
    w = hann_window(length(x))
    x=x.*w;
    if(nfft>length(x))
        x=[x;zeros(nfft-length(x))]
    end
    y = fft(x)
    #y = mag(y)
    return y
end
function fftshift(x)
    return [x[end÷2+1:end];x[1:end÷2]]
end
function read_mat(filename)
    file = matopen(filename)
    y=read(file,"samples")
    println(typeof(y))
    close(file)
    y=floor.(Int,y);
    return y
end
function submed!(x)
    x=x.-median(x)
end
function cfar!(x)
    x[:,:]=mapslices(submed!,x,dims=1)
    x[:,:]=mapslices(submed!,x,dims=2)
end
function mti!(S,P)
    # moving target indication
    # S is current frame, P is previous frame
    S_c = copy(S) # copy current frame
    S[:] = S-P # remove old frame from new
    P[:]=S_c # save copied current frame to previous frame
end
function ranging!(buffer,x,nfft,maxRange,saveStart,P)
    S=(windowed_fft(x,nfft))
    S_trim = S[saveStart:saveStart+maxRange-1]
    mti!(S_trim, P)
    pushfirst!(buffer,mag.(S_trim))
end
function doppler!(buffer,samples,nfft,maxSpeed)
    S=windowed_fft(filterSamplesDoppler(samples),nfft) 
    pushfirst!(buffer,mag.(S[1:maxSpeed]))
end
function show_waterfall(buffer,apply_cfar=false::Bool)
    p=hcat(buffer...)'
    if apply_cfar == true
        cfar!(p)
    end
    heatmap(p,c=:darkrainbow,clims=[0.,100.])
    gui()
end
function filterSamplesDoppler(dopplerSamples)
    scale_factor = 3.3/2^12
    dopplerSamples = dopplerSamples .- mean(dopplerSamples)
    dopplerSamples = dopplerSamples*scale_factor.-(3.3/2)
    return dopplerSamples
end
function preload_buffer(bufferSize,maxSpeed)
    b = CircularBuffer{Array{Float64}}(bufferSize)
    zArray = zeros(maxSpeed)
    for i=1:bufferSize
        pushfirst!(b,zArray)
    end
    return b
end
function save_mat_file(x,filename::String="sam.mat")
    # save to disk
    println("saving to ",filename)
    x=x';
    x=convert(Array{Float64,2},x)
    println("There were ",length(x)," samples")
    println(typeof(x))
    file = matopen(filename,"w")
    write(file,"samples",x)
    close(file)
    
end


#end # end module

function main(args)
    if length(args) == 2
        println("Calling ",args[1]," on ", args[2])
        getfield(liveRadar,Symbol(args[1]))(args[2])
    elseif length(args) == 3
        println("Calling ",args[1]," on port ", args[2], " at baud ", args[3])
        getfield(liveRadar,Symbol(args[1]))(args[2],parse(Int,args[3]))
    else
        println("Gotta call a function")
    end
    println("Press Enter to close")
    readline()
end

#main(ARGS)