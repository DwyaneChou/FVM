	  �  5   k820309    V          13.1        vJ`                                                                                                           
       C:\Users\Choull\Desktop\FVM\FVM\qr_solver_mod.F90 QR_SOLVER_MOD                                                    
                                                                                                           4                                                                                                    8#         @                                                      #QR_SOLVER%MATMUL    #QR_SOLVER%TRANSPOSE    #M    #N    #A 	   #B 
   #X    #X0                                                    MATMUL                                                 TRANSPOSE            @                                                    D @                                                    @                              	                    
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                                               
                    
     p          5 � p        r        5 � p        r                               D @                                                  
     p          5 � p        r        5 � p        r                                                                                   
     p          5 � p        r        5 � p        r                      #         @                                                     #GRAM_DEC%DOT_PRODUCT    #GRAM_DEC%SQRT    #A    #Q    #R    #M    #N                                                    DOT_PRODUCT                                                 SQRT          
                                                     
      p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               D                                                    
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               D                                                    
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                
                                                      
                                            #         @                                                     #UPTRI%DOT_PRODUCT    #A    #B    #X    #N                                                    DOT_PRODUCT                                                              
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                                                                                   
     p          5 � p        r        5 � p        r                               D                                                    
     p          5 � p        r        5 � p        r                                                                            #         @                                                      #CHOL_EQ%TRANSPOSE    #A    #B    #X     #N                                                    TRANSPOSE          D @                                                  
       p        5 � p        r    p          5 � p        r      5 � p        r        5 � p        r      5 � p        r                               D @                                                  
     p          5 � p        r        5 � p        r                               D @                                                   
     p          5 � p        r        5 � p        r                                D @                                          #         @                                  !                   #CHOL%DOT_PRODUCT "   #CHOL%SQRT #   #A $   #L &   #N %                                              "     DOT_PRODUCT                                            #     SQRT                                          $                    
       p        5 � p        r %   p          5 � p        r %     5 � p        r %       5 � p        r %     5 � p        r %                              D                                &                    
       p        5 � p        r %   p          5 � p        r %     5 � p        r %       5 � p        r %     5 � p        r %                                                               %            #         @                                  '                   #DOWNTRI%DOT_PRODUCT (   #A )   #B +   #X ,   #N *                                              (     DOT_PRODUCT                                          )                    
       p        5 � p        r *   p          5 � p        r *     5 � p        r *       5 � p        r *     5 � p        r *                                                              +                    
     p          5 � p        r *       5 � p        r *                              D                                ,                    
      p          5 � p        r *       5 � p        r *                                                               *            #         @                                   -                   #CG%DOT_PRODUCT .   #CG%MATMUL /   #A 0   #B 2   #X 3   #X0 4   #N 1                                              .     DOT_PRODUCT                                            /     MATMUL                                          0                    
 !      p        5 � p        r 1   p          5 � p        r 1     5 � p        r 1       5 � p        r 1     5 � p        r 1                                                              2                    
 "    p          5 � p        r 1       5 � p        r 1                              D                                3                    
 #    p          5 � p        r 1       5 � p        r 1                                                              4                    
 $    p          5 � p        r 1       5 � p        r 1                                                               1               �   H      fn#fn    �   @   J   CONSTANTS_MOD %   (  q       I_KIND+CONSTANTS_MOD %   �  q       R_KIND+CONSTANTS_MOD    
  �       QR_SOLVER !   �  ?      QR_SOLVER%MATMUL $   �  B      QR_SOLVER%TRANSPOSE    -  @   a   QR_SOLVER%M    m  @   a   QR_SOLVER%N    �  $  a   QR_SOLVER%A    �  �   a   QR_SOLVER%B    �  �   a   QR_SOLVER%X    9  �   a   QR_SOLVER%X0    �  �       GRAM_DEC %   �  D      GRAM_DEC%DOT_PRODUCT    �  =      GRAM_DEC%SQRT      $  a   GRAM_DEC%A    *	  $  a   GRAM_DEC%Q    N
  $  a   GRAM_DEC%R    r  @   a   GRAM_DEC%M    �  @   a   GRAM_DEC%N    �  {       UPTRI "   m  D      UPTRI%DOT_PRODUCT    �  $  a   UPTRI%A    �  �   a   UPTRI%B    �  �   a   UPTRI%X    =  @   a   UPTRI%N    }  {       CHOL_EQ "   �  B      CHOL_EQ%TRANSPOSE    :  $  a   CHOL_EQ%A    ^  �   a   CHOL_EQ%B      �   a   CHOL_EQ%X    �  @   a   CHOL_EQ%N      �       CHOL !   �  D      CHOL%DOT_PRODUCT    �  =      CHOL%SQRT    	  $  a   CHOL%A    -  $  a   CHOL%L    Q  @   a   CHOL%N    �  }       DOWNTRI $     D      DOWNTRI%DOT_PRODUCT    R  $  a   DOWNTRI%A    v  �   a   DOWNTRI%B    *  �   a   DOWNTRI%X    �  @   a   DOWNTRI%N      �       CG    �  D      CG%DOT_PRODUCT    �  ?      CG%MATMUL    0  $  a   CG%A    T  �   a   CG%B      �   a   CG%X    �  �   a   CG%X0    p  @   a   CG%N 