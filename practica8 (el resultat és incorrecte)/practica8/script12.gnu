# Format i nom de la imatge
set term pngcairo enhanced font 'Verdana,9'
set output "P8-1920-fig1.png"

# Permet escriure lletres gregues i altres mogudes
set encoding utf8

# Mostra els eixos
#set xzeroaxis
#set yzeroaxis

# Títol del gràfic
set title "{/Symbol F}(x)"

# Rang dels eixos
set xrange[-4:0.57]
#set yrange[-1.00e-10:1.00e10]

# Títols dels eixos
set xlabel "Distància, x (m⁻¹⁰)"
set ylabel "Funció d'ona, {/Symbol F} (m⁻⁵)"

# Canvia els nombres dels eixos per nombres personalitzats
#set ytics("1x10^-^1^0" 1.00e-10,"1x10^-^0^5" 1.00e-05,"1x10^0" 1.00e+00,"1x10^5" 1.00e+05,"1x10^1^0" 1.00e+10)
#set xtics("1x10^-^3" 1.00e-03,"1x10^-^2" 1.00e-02,"1x10^-^1" 1.00e-01,"1x10^0" 1.00e+00,"1x10^1" 1.00e+01)

# Format dels nombres dels eixos
set format y '%.2f'
set format x '%.2f'

# Escala dels eixos logarítmica
#set logscale y
#set logscale x

# Posició de la llegenda
set key top left

# Plot 
plot "aux.dat" index 0 using 1:2 with points t "{/Symbol F}(x) E1=-21eV", \
"aux.dat" index 1 using 1:2 with points t "{/Symbol F}(x) E2=-20.5eV", \
"aux.dat" index 2 using 1:2 with points t "{/Symbol F}(x) E3=-14eV", \
"aux.dat" index 3 using 1:2 with points t "{/Symbol F}(x) E4=-13eV"
#pause -1
