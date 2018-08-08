program MuscleSpindleExample
    use MuscleSpindleClass
    use MotorUnitPoolClass
    use NeuralTractClass
    use ConfigurationClass
    use ogpf 
    use randomSeedInitialize
    use MuscleNoHillClass
    use CharacterArrayClass
    use CharacterMatrixClass
    use QueueClass
    use SynapsesFactoryModule
    implicit none 
    type(Configuration) :: conf
    real(wp), parameter :: pi = 4 * atan(1.0_wp)    
    real(wp) :: dt
    real(wp) :: tf
    integer :: timeLength
    integer :: i, j
    real(wp), dimension(:), allocatable :: t, MNv_mV, l, v, a
    real(wp) :: tic, toc
    type(gpf) :: gp
    character(len = 80) :: muscle = 'TA'
    character(len = 80) :: filename = 'confMuscleSpindle.rmto'
    type(CharacterMatrix) :: Synapses
    type(MuscleSpindle):: spindle
    real(wp) :: Amp, f
    real(wp), dimension(:), allocatable :: IaFR, fBag1, fBag2, fChain
    real(wp), dimension(:), allocatable :: TBag1, TBag2, TChain

    conf = Configuration(filename)

    spindle = MuscleSpindle(conf, muscle)
    
    dt = conf%timeStep_ms
    tf = conf%simDuration_ms
    timeLength = int(tf/dt)

    t = [(dt*(i-1), i=1, timeLength)]


    allocate(l(timeLength))
    allocate(v(timeLength))
    allocate(a(timeLength))
    allocate(IaFR(timeLength))
    allocate(fBag1(timeLength))
    allocate(fBag2(timeLength))
    allocate(fChain(timeLength))
    allocate(TBag1(timeLength))
    allocate(TBag2(timeLength))
    allocate(TChain(timeLength))



    Amp = 0.005
    f = 0.25


    l = 0.995 - Amp*cos(2*pi*f*t/1000.0)
    v = 2*pi*f*Amp*sin(2*pi*f*t/1000.0)/1000.0
    a = 4*(pi**2)*(f**2)*Amp*cos(2*pi*f*t/1000.0)/(1000.0**2)
    
    call cpu_time(tic)
    do i = 1, timeLength
        call spindle%atualizeMuscleSpindle(t(i), l(i), v(i), a(i), 10.0_wp,10.0_wp)
        IaFR(i) = spindle%IaFR_Hz
        fBag1(i) = spindle%fusimotorActivation(1)
        fBag2(i) = spindle%fusimotorActivation(2)
        fChain(i) = spindle%fusimotorActivation(3)
        TBag1(i) = spindle%fiberTension(1)
        TBag2(i) = spindle%fiberTension(3)
        TChain(i) = spindle%fiberTension(5)
    end do
    call cpu_time(toc)

    print '(F15.6, A)', toc - tic, ' seconds'

    call gp%title('Bag1 Firing Rate')
    call gp%xlabel('t (ms))')
    call gp%ylabel('FR Bag1 (Hz)')
    call gp%plot(t, fBag1, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Ia Firing Rate')
    call gp%xlabel('t (ms))')
    call gp%ylabel('Ia FR (Hz)')
    call gp%plot(t, IaFR, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Ia Firing Rate')
    call gp%xlabel('fascile Length (m)')
    call gp%ylabel('Ia FR (Hz)')
    call gp%plot(l, IaFR, 'with line lw 2 lc rgb "#0008B0"')

    call gp%title('Ia Firing Rate')
    call gp%xlabel('Fascicle velocity (m/ms)')
    call gp%ylabel('Ia FR (Hz)')
    call gp%plot(v, IaFR, 'with line lw 2 lc rgb "#0008B0"')
end program MuscleSpindleExample