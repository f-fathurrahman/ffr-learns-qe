module nvtx_m
use iso_c_binding
implicit none
integer,parameter :: num = 20
integer,private :: col(num) = &
    [int(Z'0000ff00'), &
     int(Z'00ffe4c4'), &
     int(Z'00ffa500'), &
     int(Z'0000ffff'), &
     int(Z'00ee82ee'), &
     int(Z'00ffff00'), &
     int(Z'00006400'), &
     int(Z'00808000'), &
     int(Z'00ee00ee'), &
     int(Z'00008080'), &
     int(Z'00ffd700'), &
     int(Z'006a5acd'), &
     int(Z'00dda0dd'), &
     int(Z'00ffc0cb'), &
     int(Z'00808080'), &
     int(Z'000000ff'), &
     int(Z'00d2691e'), &
     int(Z'00ff0000'), &
     int(Z'00cd853f'), &
     int(Z'00ccccff')]
character(len=256),private :: tempName
type, bind(C):: nvtxEventAttributes
  integer(C_INT16_T):: version=1
  integer(C_INT16_T):: size=48 !
  integer(C_INT):: category=0
  integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
  integer(C_INT):: color
  integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
  integer(C_INT):: reserved0
  integer(C_INT64_T):: payload ! union uint,int,double
  integer(C_INT):: messageType=1 ! NVTX_MESSAGE_TYPE_ASCII = 1
  type(C_PTR):: message ! ascii char
end type
contains
subroutine nvtxStartRange(name,id)
  character(kind=c_char,len=*) :: name
  integer, optional:: id
end subroutine
subroutine nvtxEndRange
end subroutine
end module nvtx_m
