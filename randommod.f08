  ! randommod.f08 :: Module for pseudo random number generation. The internal pseudo random
  !                  generator is the xoroshiro128plus method.
  !
  MODULE randommod

    implicit none
  
    private

      ! A 64 bit floating point type
      integer, parameter :: dp = KIND(0.0d0)

      ! A 32 bit integer type
      integer, parameter :: i4 = SELECTED_INT_KIND(9)

      ! A 64 bit integer type
      integer, parameter :: i8 = SELECTED_INT_KIND(18)

    public
      
      !
      ! DESCRIPTION :: Random number generator type, which contains the state
      !
      TYPE rng_t
  
          ! The rng state (always use your own seed)
          integer(i8), private :: s(2) = [123456789_i8, 987654321_i8]
          integer(i8), private :: separator(32) ! Separate cache lines (parallel use)
   
        CONTAINS
   
          procedure, non_overridable :: set_seed    ! Seed the generator
          procedure, non_overridable :: jump        ! Jump function (see below)
          procedure, non_overridable :: int_4       ! 4-byte random integer
          procedure, non_overridable :: int_8       ! 8-byte random integer
          procedure, non_overridable :: unif_01     ! Uniform (0,1] real
          procedure, non_overridable :: two_normals ! Two normal(0,1) samples
          procedure, non_overridable :: poisson     ! Sample from Poisson-dist.
          procedure, non_overridable :: circle      ! Sample on a circle
          procedure, non_overridable :: sphere      ! Sample on a sphere
          procedure, non_overridable :: next        ! Internal method
     
      END TYPE rng_t

      !
      ! DESCRIPTION :: Parallel random number generator type
      !
      TYPE prng_t
     
        type(rng_t), allocatable :: rngs(:)
      
        CONTAINS
   
          procedure, non_overridable :: init_parallel
          procedure, non_overridable :: update_seed
      
      END TYPE prng_t

    CONTAINS

      !
      ! DESCRIPTION :: Initialize a collection of rng's for parallel use
      ! PARAMETERS :: class(prng_t) self, integer n_proc, class(rng_t) rng
      !
      SUBROUTINE init_parallel(self, n_proc, rng)
    
        class(prng_t), intent(inout) :: self
        type(rng_t), intent(inout)   :: rng
        integer, intent(in)          :: n_proc
        integer                      :: n

        if (n_proc < 1) error stop "init_parallel: n_proc < 1"

        allocate(self%rngs(n_proc))
        self%rngs(1) = rng

        do n = 2, n_proc
          self%rngs(n) = self%rngs(n-1)
          call self%rngs(n)%jump()
        end do
    
      END SUBROUTINE init_parallel

      !
      ! DESCRIPTION :: Parallel RNG instances are often used temporarily. This routine can
      !                afterwards be used to update the seed of the user's sequential RNG.
      ! PARAMETERS  :: class(prng_t) self, class(rng_t) rng
      !
      SUBROUTINE update_seed(self, rng)
    
        class(prng_t), intent(inout) :: self
        type(rng_t), intent(inout)   :: rng
        integer                      :: n

        do n = 1, size(self%rngs)
          ! Perform exclusive-or with each parallel rng
          rng%s(1) = ieor(rng%s(1), self%rngs(n)%s(1))
          rng%s(2) = ieor(rng%s(2), self%rngs(n)%s(2))
        end do
    
      END SUBROUTINE update_seed

      !
      ! DESCRIPTION :: Set a seed for the rng
      ! PARAMETERS  :: class(prng_t) self, integer the_seed
      !
      SUBROUTINE set_seed(self, the_seed)
    
        class(rng_t), intent(inout) :: self
        integer(i8), intent(in)     :: the_seed(2)

        self%s = the_seed

        ! Simulate calls to next() to improve randomness of first number
        call self%jump()
    
      END SUBROUTINE set_seed

      !
      ! DESCRIPTION :: This is the jump function for the generator. It is equivalent
      !                to 2^64 calls to next(); it can be used to generate 2^64
      !                non-overlapping subsequences for parallel computations.
      ! PARAMETERS  :: class(prng_t) self
      !
      SUBROUTINE jump(self)
  
        class(rng_t), intent(inout) :: self
        integer                     :: i, b
        integer(i8)                 :: t(2), dummy

        ! The signed equivalent of the unsigned constants
        integer(i8), parameter :: jmp_c(2) = (/-4707382666127344949_i8, -2852180941702784734_i8/)

        t = 0
        do i = 1, 2
          do b = 0, 63
            if (iand(jmp_c(i), shiftl(1_i8, b)) /= 0) then
              t = ieor(t, self%s)
            end if
          
            dummy = self%next()
        
          end do
        end do

        self%s = t
    
      END SUBROUTINE jump

      !
      ! DESCRIPTION :: Return 4-byte integer
      ! PARAMETERS  :: class(rng_t) self
      ! RETURNS     :: integer(i4) y
      !
      FUNCTION int_4(self) RESULT(y)
  
        class(rng_t), intent(inout) :: self
        integer(i4) :: y
    
        y = int( self%next(), i4 )
    
      END FUNCTION int_4

      !
      ! DESCRIPTION :: Return 8-byte integer
      ! PARAMETERS  :: class(rng_t) self 
      ! RETURNS     :: integer(i8) y 
      !
      FUNCTION int_8(self) RESULT(y)
        
        class(rng_t), intent(inout) :: self
        integer(i8) :: y
        
        y = self%next()
      
      END FUNCTION int_8

      !
      ! DESCRIPTION :: Get a uniform [0,1) random real (double precision)
      ! PARAMETERS  :: class(rng_t) self
      ! RETURNS     :: real(dp) y
      !
      FUNCTION unif_01(self) RESULT(y)
        
        class(rng_t), intent(inout) :: self
        integer(i8)                 :: x
        real(dp)                    :: tmp, y
        
        x = self%next()
        x = ior(shiftl(1023_i8, 52), shiftr(x, 12))
        
        y = transfer(x, tmp) - 1.0_dp
      
      END FUNCTION unif_01

      !
      ! DESCRIPTION  :: Return two normal random variates with mean 0 and variance 1.
      !                 http://en.wikipedia.org/wiki/Marsaglia_polar_method
      ! PARAMETERS   :: class(rng_t) self
      ! RETURNS      :: real(dp) y(2)
      !
      FUNCTION two_normals(self) RESULT(rands)
        
        class(rng_t), intent(inout) :: self
        real(dp)                    :: y(2), sum_sq

        do
        
          y(1) = 2 * self%unif_01() - 1
          y(2) = 2 * self%unif_01() - 1
          
          sum_sq = sum(y**2)
          
          if( sum_sq < 1.0_dp .and. sum_sq > 0.0_dp ) exit
        
        end do
        
        y = y * sqrt(-2 * log(sum_sq) / sum_sq)
      
      END FUNCTION two_normals

      !
      ! DESCRIPTION :: Return Poisson random variate with rate lambda. Works well for lambda < 30
      !                or so. For lambda >> 1 it can produce wrong results due to roundoff error.
      ! PARAMETERS  :: class(rng_t) self, real(dp) lambda
      ! RETURNS     :: integer(i4) y
      !
      FUNCTION poisson(self, lambda) RESULT(y)
        
        class(rng_t), intent(inout) :: self
        real(dp), intent(in)        :: lambda
        integer(i4)                 :: y
        real(dp)                    :: expl, p

        expl = exp(-lambda)
        y    = 0
        p    = self%unif_01()

        do while( p > expl )
          y = y + 1
          p = p * self%unif_01()
        end do
      
      END FUNCTION poisson

      !
      ! DESCRIPTION :: Sample point on a circle with given radius
      ! PARAMETERS  :: class(rng_t) self
      ! RETURNS     :: real(dp) y(2)
      !
      FUNCTION circle(self, radius) RESULT(y)
        
        class(rng_t), intent(inout) :: self
        real(dp), intent(in)        :: radius
        real(dp)                    :: rands(2), y(2)
        real(dp)                    :: sum_sq

        ! Method for uniform sampling on circle
        do
        
          rands(1) = 2 * self%unif_01() - 1
          rands(2) = 2 * self%unif_01() - 1
          sum_sq   = sum(rands**2)
          if (sum_sq <= 1) exit
        
        end do

        y(1) = (rands(1)**2 - rands(2)**2) / sum_sq
        y(2) = 2 * rands(1) * rands(2) / sum_sq
        y    = y * radius
      
      END FUNCTION circle

      !
      ! DESCRIPTION :: Sample point on a sphere with given radius
      ! PARAMETERS  :: class(rng_t) self, real(dp) radius
      ! RETURNS     :: real(dp) y(3)
      !
      function sphere(self, radius) result(y)
        
        class(rng_t), intent(inout) :: self
        real(dp), intent(in)        :: radius
        real(dp)                    :: rands(2), y(3)
        real(dp)                    :: sum_sq, tmp_sqrt

        ! Marsaglia method for uniform sampling on sphere
        do
        
          rands(1) = 2 * self%unif_01() - 1
          rands(2) = 2 * self%unif_01() - 1
          sum_sq   = sum(rands**2)
          if (sum_sq <= 1) exit
        
        end do

        tmp_sqrt = sqrt(1 - sum_sq)
        y(1:2) = 2 * rands(1:2) * tmp_sqrt
        y(3)   = 1 - 2 * sum_sq
        y      = y * radius
      
      END FUNCTION sphere

      !
      ! DESCRIPTION :: Interal routine - get the next value (returned as 64 bit signed integer)
      ! PARAMETERS  :: class(rng_t) self
      ! RETURNS     :: integer(i8) y
      !
      FUNCTION next(self) RESULT(y)
        
        class(rng_t), intent(inout) :: self
        integer(i8)                 :: y
        integer(i8)                 :: t(2)

        t         = self%s
        y         = t(1) + t(2)
        t(2)      = ieor(t(1), t(2))
        self%s(1) = ieor(ieor(rotl(t(1), 55), t(2)), shiftl(t(2), 14))
        self%s(2) = rotl(t(2), 36)
      
      END FUNCTION next

      !
      ! DESCRIPTION :: Helper function for next()
      ! PARAMETERS  :: integer(8) x, integer k
      ! RETURNS     :: integer(i8) y
      !
      PURE FUNCTION rotl(x, k) RESULT(y)
        
        integer(i8), intent(in) :: x
        integer, intent(in)     :: k
        integer(i8)             :: y

        y = ior(shiftl(x, k), shiftr(x, 64 - k))
      
      END FUNCTION rotl

  END MODULE randommod
