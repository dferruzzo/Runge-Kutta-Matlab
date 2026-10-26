# Runge-Kutta em MATLAB

Implementacoes de metodos de Runge-Kutta para integracao numerica de
problemas de valor inicial, acompanhadas por exemplos com o atrator de Lorenz.

## Funcoes

### `rk4`

Executa um passo do metodo classico de Runge-Kutta de quarta ordem.

```matlab
[t, y] = rk4(f, t0, y0, h)
```

- `f`: funcao da EDO no formato `f(t, y)`;
- `t0`, `y0`: tempo e estado iniciais;
- `h`: tamanho do passo;
- `t`, `y`: tempo e estado apos um passo.

### `rk45`

Executa um passo do metodo de Runge-Kutta-Fehlberg. O estado retornado usa a
aproximacao de quarta ordem e a funcao tambem devolve a estimativa de erro
obtida pelo par incorporado de quarta e quinta ordem.

```matlab
[t, y, error] = rk45(f, t0, y0, h)
```

- `error`: estimativa vetorial do erro local do passo.

### `rk45_h_adj`

Executa um passo de Runge-Kutta-Fehlberg e calcula um proximo passo sugerido
a partir de uma tolerancia.

```matlab
[t, y, error, h_adj] = rk45_h_adj(f, t0, y0, h, epsilon)
```

- `epsilon`: tolerancia usada no ajuste;
- `h_adj`: tamanho de passo sugerido para a proxima iteracao.

O chamador deve limitar `h_adj` aos valores minimo e maximo desejados e
rejeitar passos quando necessario.

## Exemplos

- `tests.m`: compara `rk4` e `rk45` na simulacao do atrator de Lorenz.
- `tests2.m`: integra o atrator de Lorenz com `rk45_h_adj`, limitando o passo
  entre `h_min` e `h_max`.

No MATLAB ou GNU Octave, execute um dos scripts a partir da raiz do projeto:

```matlab
tests
% ou
tests2
```

## Referencia

Ramin Esfandiari, *Numerical Methods for Engineers and Scientists Using
MATLAB*, 2013.
