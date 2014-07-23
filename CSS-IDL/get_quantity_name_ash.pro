function get_quantity_name, quantity, tex=tex

If (keyword_set(tex)) Then Begin
   case quantity of

     1    : name=textoidl('V_r')
     2    : name=textoidl('V_{\theta}')
     3    : name=textoidl('V_{\phi}')
     4    : name='Entropy'
     5    : name='Pressure'
     6    : name="T'"
     7    : name=textoidl('rho u_r')
     8    : name=textoidl('rho u_{\theta}')
     9    : name=textoidl('rho u_{\phi}')
     10   : name='Kinetic Energy'
     11   : name=textoidl('\theta Kinetic Energy Flux')
     12   : name='Radial Kinetic Energy Flux'
     13   : name=textoidl('Momentum Flux: r \theta')
     14   : name=textoidl('Momentum Flux: r \phi')
     15   : name=textoidl('Momentum Flux: \theta \phi')
     16   : name=textoidl('\theta Enthalpy Flux')
     17   : name='Radial Enthalpy Flux'
     18   : name='Pressure Work'
     19   : name='Kinetic Helicity'
     20   : name='Enstrophy'
     21   : name='Bouancy Work'
     22   : name='Density'
     24   : name='Entropy Gradient'
     25   : name='Ra'
     26   : name='DRKE'
     27   : name='MCKE'
     28   : name=textoidl('F_{un}' )
     29   : name=textoidl('F_{\rad}')
     30   : name='Ro'
     31   : name='Conv. Ro'
     32   : name=textoidl('\theta Viscous Flux')
     33   : name='Radial Viscous Flux'
     34   : name=textoidl('\theta Potential Energy Flux')
     35   : name='Radial Potential Energy Flux'
     36   : name='Viscous Work'
     37   : name='NAKE'
     43   : name=textoidl('B_r')
     44   : name=textoidl('B_{\theta}')
     45   : name=textoidl('B_{\phi}')
     46   : name=textoidl('J_r')
     47   : name=textoidl('J_{\theta}')
     48   : name=textoidl('J_{\phi}')
     49   : name='Magnetic Energy'
     50   : name=textotidl('\theta Mag. Energy Flux')
     51   : name='Radial Mag. Energy Flux'
     58   : name='Magnetic Helicity'
     59   : name='Rate of Magnetic Work'
     60   : name='Cross Helicity'
     61   : name=textoidl('J^2')
     62   : name='Current Helicity'
     63   : name='TME'
     64   : name='PME'
     65   : name='NAME'
     66   : name='NATME'
     67   : name='NAPME'
     else : name=' '

   endcase
Endif Else Begin
   case quantity of

     1    : name='Vr'
     2    : name='Vtheta'
     3    : name='Vphi'
     4    : name='Entropy'
     5    : name='Pressure'
     6    : name='Temperature'
     7    : name='Radial Mass Flux'
     8    : name='Theta Mass Flux'
     9    : name='Phi Mass Flux'
     10   : name='Kinetic Energy'
     11   : name='Theta Kinetic Energy Flux'
     12   : name='Radial Kinetic Energy Flux'
     13   : name='Momentum Flux: r theta'
     14   : name='Momentum Flux: r phi'
     15   : name='Momentum Flux: theta phi'
     16   : name='Momentum Flux: r r'
     17   : name='Momentum Flux: theta theta'
     18   : name='Momentum Flux: phi phi'
     19   : name='Theta Enthalpy Flux'
     20   : name='Radial Enthalpy Flux'
     21   : name='Pressure Work'
     22   : name='Kinetic Helicity'
     23   : name='Enstrophy'
     24   : name='Bouancy Work'
     25   : name='Radial Diffusive Flux'
     26   : name='Theta Diffusive Flux'
     27   : name='Phi Disffusive Flux'
     28   : name='Density'
     29   : name='Viscous Work'
     30   : name='KE Advection'
     31   : name='Theta Viscous Flux'
     32   : name='Radial Viscous Flux'
     33   : name='Theta Potential Energy Flux'
     34   : name='Radial Potential Energy Flux'
     35   : name='Theta Entropy Flux'
     36   : name='Radial Entropy Flux'
     37   : name='Theta Radiative Flux'
     38   : name='Radial Radiative Flux'
     39   : name='Total Flux'
     40   : name='Angular Momentum'
     41   : name='ds/dr'
     48   : name='s_pert'
     49   : name='ds/dr_pert'
     50   : name='omega_r'
     51   : name='omega_theta'
     52   : name='omega_phi'
     61   : name='Br'
     62   : name='Btheta'
     63   : name='Bphi'
     64   : name='Radial Current'
     65   : name='Theta Current'
     66   : name='Phi Current'
     67   : name='Magnetic Energy'
     68   : name='Theta Poynting Flux'
     69   : name='Radial Poynting Flux'
     70   : name='mixed flux (ut br)'
     71   : name='mixed flux (ut bt)'
     72   : name='mixed flux (ut bp)'
     73   : name='mixed flux (ur br)'
     74   : name='mixed flux (ur bt)'
     75   : name='mixed flux (ur bp)'
     76   : name='Magnetic Helicity'
     77   : name='Magnetic Work'
     78   : name='Cross Helicity'
     79   : name='J^2'
     80   : name='Current Helicity'
     81   : name='Omega effect'
     82   : name='Alpha effect'
     83   : name='radial mean emf'
     84   : name='theta mean emf'
     85   : name='phi mean emf'
     86   : name='mean Bp prod. term.'
     87   : name='fluc. Bp prod. term.'
     110  : name='dynamo shear r'
     111  : name='dynamo shear t'
     112  : name='dynamo shear p'
     113  : name='mean dynamo shear r'
     114  : name='mean dynamo shear t'
     115  : name='mean dynamo shear p'
     116  : name='radial dynamo adv.'
     117  : name='theta dynamo adv.'
     118  : name='phi dynamo adv.'
     119  : name='mean radial dynamo adv.'
     120  : name='mean theta dynamo adv.'
     121  : name='mean phi dynamo adv.'
     122  : name='radial dynamo comp. term'
     123  : name='theta dynamo comp. term'
     124  : name='phi dynamo comp. term'
     125  : name='radial dynamo diff. term'
     126  : name='theta dynamo diff. term'
     127  : name='phi dynamo diff. term'
     128  : name='mag SGS r'
     129  : name='mag SGS t'
     130  : name='mag SGS p'
     131  : name='mag. (br bt)'
     132  : name='mag. (br bp)'
     133  : name='mag. (bt bp)'
     134  : name='mag. (br br)'
     135  : name='mag. (bt bt)'
     136  : name='mag. (bp bp)'
     201  : name='EFB_visc_r'
     202  : name='EFB_visc_t'
     203  : name='EFB_sdif_r'
     204  : name='EFB_sdif_t'
     205  : name='EFB_tdif_r'
     206  : name='EFB_tdif_t'
     207  : name='EFB_ketb_r'
     208  : name='EFB_ketb_t'
     209  : name='EFB_kemn_r'
     210  : name='EFB_kemn_t'
     211  : name='EFB_entb_r'
     212  : name='EFB_entb_t'
     213  : name='EFB_enmn_r'
     214  : name='EFB_enmn_t'

     else : name='Unknown'

   endcase
EndElse
   return, name

end






