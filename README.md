# Tarea nro. 7
## FI3104B - Métodos Numéricos para la Ciencia y la Ingeniería
#### Prof. Valentino González


Integre las ecuaciones para &rho; y _v_ en un fluido compresible en una
dimension, sin viscosidad ni gravedad, y usando el método de las características
como fue visto en clase. Las ecuaciones a resolver son:

<img src="eqs/fluido.png" height="110px"/>

> Latex `\begin{flalign*} \dfrac{\partial \rho}{\partial t} + \dfrac{\partial \rho v}{\partial x} &= 0\\ \rho \dfrac{\partial v}{\partial t} + \rho v \dfrac{\partial v}{\partial x} &= - \dfrac{\partial P}{\partial x}\\ P &= A \rho^\gamma \end{flalign*}`


La última es la ecuación de estado. Para ella utilice &gamma; = 5/3 (a qué
corresponde?) con A(x &le; 1/3) = 4 y A(x > 1/3) = 1.

> NOTA: el cambio de ecuación de estado en x = 1/3 hará que los pulsos que se
> propagan por el fluido se reflejen y difracten al pasar por la discontinuidad.


Considere x entre 0 y 1 con las siguientes condiciones iniciales y de borde:

<img src="eqs/iniciales.png" height="110px"/>

> Latex `\begin{flalign*} v(0, t) &= 0\\ v(1, t) &= 0\\ \rho(x > 0.1, 0) &= 1\\ \rho(x \leq 0.1, 0) & = 1 + 0.0RRR (1 + cos(10 \pi x))\\ \end{flalign*}`

> NOTA: la condición inicial es un pulso de sobredensidad en x=0 que se
> propagará por el fluido.

donde RRR son los últimos 3 dígitos de su RUT (antes del dígito verificador).
Integre las ecuaciones hasta que el tiempo alcance al menos un valor de 1.2 al
centro del intervalo 0 &leq; x &leq; 1.

Para resolver el problema divida el intervalo 0 &leq; x &leq; 1 en a lo menos 10
mil trazos iguales. Debe diseñar con cuidado el algoritmo para tratar los
bordes. Explique bien en su informe lo que hizo. También recuerde que el
algoritmo es levemente distinto al pasar de tiempo par a tiempo impar y
vice-versa.

> NOTA: A diferencia de tareas anteriores, el output de su código será mucho más
> masivo, puede que tarde más en correr y no tendrá necesariamente una grilla
> uniforme en x y t. Por estos motivos se recomienda que guarde el output en un
> archivo separado y escriba un script aparte para crear los gráficos. Más
> detalles a continuación.

1. Dibuje las funciones &rho; y _v_ como función de x para t = 0, 0.2, 0.4, 0.6,
   0.8, y 1.0. Para cada uno de esos tiempo haga un solo plot que combine
   &rho;(x) y (1 + _v_(x)).

1. En el plano (x, t), dibuje las superficie &rho;(x, t) y _v_(x, t) de modo que
   se pueda ver la reflexión que ocurre en x = 1/3. Esto lo puede hacer usando
   curvas de nivel o ploteos de superficie como hizo en las tareas anteriores.

> NOTAS PARA LOS GRAFICOS
>
> En su algoritmo Ud. implementará un barrido sobre x luego del cual se habrá
> avanzado en un paso temporal. Luego de cada paso temporal tendrá el mismo
> número de valores x y t, pero los pasos ya no serán uniformes ni en x ni en t.
> La recomendación es que no guarde todos los valores, sino que cada, por
> ejemplo, 100 barridos, guarde 1 de cada 100 valores de x, t, &rho; y _v_.
> Luego utilice rutinas como `plt.pcolormesh` para hacer los gráficos de
> superficie.
>
> Para los otros gráficos, aquellos a tiempo constante, se sugiere usar alguna
> rutina de interpolación como `scipy.interpolate.interp2d` para obtener los
> valores de &rho; y _v_ en un tiempo dado, para todos los x.

__Otras Notas.__

- Utilice `git` durante el desarrollo de la tarea para mantener un historial de
  los cambios realizados. La siguiente [*cheat
  sheet*](https://education.github.com/git-cheat-sheet-education.pdf) le puede
  ser útil. Evaluar el uso efectivo de `git`. Recuerde hacer cambios
  significativos pero relativamente pequeños y guardar seguido.  Evite hacer
  `commits` de código que no compila y deje mensajes que permitan entender los
  cambios realizados.

- Evaluaremos su uso correcto de python. Si define una función relativamente
  larga o con muchos parámetros, recuerde escribir el *doctsring* que describa
  los parametros y que es lo que hace la función.  También recuerde usar nombres
  explicativos para las variables y las funciones.  El mejor nombre es aquel que
  permite entender que hace la función sin tener que leer su implementación.

- Los códigos entregados deben pasar las reglas de
  [PEP8](https://www.python.org/dev/peps/pep-0008/). En la línea de comando
  puede intentar `pep8 <script.py>` y asegurarse de que no hay errores ni
  advertencias. Los comandos `pep8 --show-source <script.py>` y `pep8
  --show-pep8 <script.py>` que vimos en clase le pueden ser útiles. Si es de
  aquellos que resuelven su tarea usando el `ipython notebook`, entonces exporte
  la tarea a un archivo normal de `python` y asegúrese de que pasa el test de
  PEP8.

- La tarea se entrega como un *pull request* en github. El *pull request* debe
  incluir todos los códigos usados además de su informe.

- El informe debe ser entregado en formato *pdf*, este debe ser claro sin
  información ni de más ni de menos.
