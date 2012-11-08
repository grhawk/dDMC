MODULE precision 
  
  integer, parameter :: kr4 = selected_real_kind(6,37)      
  integer, parameter :: kr8 = selected_real_kind(15,307)  
  integer, parameter :: kr16 = selected_real_kind(30,1000) 
  integer, parameter :: ki4 = selected_int_kind(9)           
  integer, parameter :: ki8 = selected_int_kind(18)          
  integer, parameter :: kc4 = kr4                            
  integer, parameter :: kc8 = kr8                            
  ! generic kinds
  integer, parameter :: ki=ki4,kr=kr8,kc=kc8,kch=64
  
END MODULE precision
