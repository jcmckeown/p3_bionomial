## Illustrating the Binomial Theorem

Draws a rotating n-dimensional simplex subdivided into k-by-(n-k)-dimensional products of simplices, which illustrates the (analytic) binomial theorem when written
$$ \frac{(a+b)^n}{n!} = \sum \frac{a^k}{k!} \frac{b^(n-k)}{(n-k)!} $$
The rotation can be guided along a 2-d foliation of SO(n) by moving the Pointing Device.

### Known omission

 * The Processing3 documentation indicates that when the Main Pointer is outside of the Sketch Window, then the sketch should see MouseX==MouseY==0;
experimentally, however, this need not be the case.  Still pondering whether to emulate said behaviour.

## Dependencies; Running/Compiling

Runs in [Processing3](https://github.com/processing/processing), requires library
[interfascia](https://github.com/brendanberg/interfascia)

P3, of course, has a built-in compiler.

## License

As far as is possible, this project is released CC-BY-SA

As far as the author is concerned, the -BY- is sufficiently satisfied by naming him in one included file.

© Jesse C. McKeown 2021
