As imagens para o relatório foram obtidas respectivamente com os comandos

python3 make_launch.py 320 90 360 -m 0.5 -d 5.5e-4 -s 15
python3 make_launch.py 891 90 360 -m 0.042 -d 2.8e-5 -t 52 -s 15

Para ajuda

python3 make_launch.py --help

Em resumo, o programa recebe como parâmetros obrigatórios a velocidade inicial,
o angulo azimutal do disparo e a direção do disparo, dada por um angulo polar.

Como parâmetros opcionais, é possível customizar: 

-l  --latitude              a latitude do disparo

-m  --mass                  a massa do projétil

-d  --drag                  o coeficiente da força de arrasto aerodinâmico

-t  --time                  o intervalo de tempo no qual o programa calculará a trajetória

-s  --scale                 float entre 0 e 100.
                            0 focaliza na trajetória com resistência do ar
                            100 focaliza na trajetória sem resistência do ar
                            é útil quando as trajetórias são muito discrepantes
                            
                            só se aplica aos gráficos em 3D e ao plot da
                            projeção da trajetória no solo                    

Em caso de dúvidas ou problemas, entre em contato.

Davi Feliciano - dfeliciano37@gmail.com

