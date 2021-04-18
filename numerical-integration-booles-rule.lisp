
;; Date: 2021-04-14
;; 
;; This program implements both Boole's rule and the composite Boole's
;; rule for integration.
;; 
;; Source:
;;   -> "https://www.scipress.com/BSMaSS.2.1.pdf"
;;       o Introduces both Boole's rule and the composite Boole's rule.
;;   -> "https://atozmath.com/example/CONM/NumeInte.aspx?he=e&q=BR"
;;       o Introduces the composite Boole's rule and also offers
;;         examples.
;;   -> "https://en.wikipedia.org/wiki/Boole%27s_rule"
;;   -> "https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas"
;;   -> "https://mathworld.wolfram.com/BoolesRule.html"
;;   -> "http://netedu.xauat.edu.cn/jpkc/netedu/jpkc2009/jsff/content/syjx/3/Chapter4.pdf"
;;   -> "https://www.whitman.edu/mathematics/calculus_late_online/section10.05.html"
;;       o Offers examples for numerical integration methods with solutions.
;;   -> "https://scicomp.stackexchange.com/questions/29701/booles-rule-in-python"
;;       o Implementation of several numerical integration methods,
;;         including Boole's rule, in the Python programming language.



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; -- Boole's rule.                                                -- ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun integrate-by-boole (f x1 x5)
  "Computes and returns the numerical integral of the function F in the
   closed interval [X1, X5], where X1 demarcates the lower bound
   of the abscissa and X5 its upper bound, across which Boole's rule is
   applied."
  (declare (type (function (real) real) f))
  (declare (type real                   x1))
  (declare (type real                   x5))
  (let ((h (/ (- x5 x1) 4)))
    (declare (type real h))
    (let* ((x2 (+ x1 h))
           (x3 (+ x2 h))
           (x4 (+ x3 h)))
      (declare (type real x2 x3 x4))
      (* (/ (* 2 h) 45)
         (+ (*  7 (funcall f x1))
            (* 32 (funcall f x2))
            (* 12 (funcall f x3))
            (* 32 (funcall f x4))
            (*  7 (funcall f x5)))))))

;;; -------------------------------------------------------

(float
  (integrate-by-boole
    #'(lambda (x) (/ 1 x))
    1 2))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; -- Boole's rule.                                                -- ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun integrate-by-boole (f a b n)
  "Computes and returns the numerical integral of the function F in the
   closed interval [A, B], divided into N equidistant subintervals, and
   obtained by applying Boole's rule to each of the N sections."
  (declare (type (function (real) real) f))
  (declare (type real                   a))
  (declare (type real                   b))
  (declare (type (integer 1 *)          n))
  ;; The function f betwixt the bounds [a, b] is divided into patches
  ;; of the width k.
  ;; Each such patch is integrated employing Boole's rule with a
  ;; distance of its points being equal to h.
  (let* ((k (/ (- b a) n))
         (h (/ k 4)))
    (declare (type real k))
    (declare (type real h))
    (loop
      repeat n
      for x1 = a then (+ x1 k)
      sum (let* ((x2 (+ x1 h))
                 (x3 (+ x2 h))
                 (x4 (+ x3 h))
                 (x5 (+ x4 h)))
            (* (/ (* 2 h) 45)
               (+ (*  7 (funcall f x1))
                  (* 32 (funcall f x2))
                  (* 12 (funcall f x3))
                  (* 32 (funcall f x4))
                  (*  7 (funcall f x5))))))))

;;; -------------------------------------------------------

(float
  (integrate-by-boole
    #'(lambda (x) (/ 1 x))
    1 2 8))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; -- Boole's rule with auxiliary function.                        -- ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun calculate-boole-area (f x1 h)
  "Computes and returns the numerical integral of the function F in the
   closed interval [X1, X1 + H], where X1 demarcates the lower bound
   of the abscissa and H the width of each division of the interval,
   across which Boole's rule is applied."
  (declare (type (function (real) real) f))
  (declare (type real                   x1))
  (declare (type real                   h))
  (let* ((x2 (+ x1 h))
         (x3 (+ x2 h))
         (x4 (+ x3 h))
         (x5 (+ x4 h)))
    (declare (type real x2 x3 x4 x5))
    (* (/ (* 2 h) 45)
       (+ (*  7 (funcall f x1))
          (* 32 (funcall f x2))
          (* 12 (funcall f x3))
          (* 32 (funcall f x4))
          (*  7 (funcall f x5))))))

(defun integrate-by-boole (f a b n)
  "Computes and returns the numerical integral of the function F in the
   closed interval [A, B], divided into N equidistant subintervals, and
   obtained by applying Boole's rule to each of the N sections."
  (declare (type (function (real) real) f))
  (declare (type real                   a))
  (declare (type real                   b))
  (declare (type (integer 1 *)          n))
  (let* ((k (/ (- b a) n))
         (h (/ k 4)))
    (declare (type real k))
    (declare (type real h))
    (loop
      repeat n
      for x1 of-type real = a then (+ x1 k)
      sum (calculate-boole-area f x1 h))))

;;; -------------------------------------------------------

(float
  (integrate-by-boole
    #'(lambda (x) (/ 1 x))
    1 2 8))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; -- Boole's rule with auxiliary function.                        -- ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; This version of the program utilizes a genuine representation of
;; Boole's rule for a single patch, here norned ``calculate-boole-area'',
;; which accepts the function to integrate and the patch's bounds. While
;; this method reflects the principle very naturally, especially in
;; terms of the arguments and its independence from the repeatingly
;; invoking ``integrate-by-boole'' function, the patch subinterval width
;; h, equal across all such patches, must be redundantly recomputed for
;; each call to ``calculate-boole-area'' by the same function, lacking
;; an argument with respect to fidelity. This imposes computational an
;; unnecessary computational penalty.

(defun calculate-boole-area (f x1 x5)
  "Computes and returns the numerical integral of the function F in the
   closed interval [X1, X5], where X1 demarcates the lower bound
   of the abscissa and X5 its upper bound, across which Boole's rule is
   applied."
  (declare (type (function (real) real) f))
  (declare (type real                   x1))
  (declare (type real                   x5))
  (let ((h (/ (- x5 x1) 4)))
    (declare (type real h))
    (let* ((x2 (+ x1 h))
           (x3 (+ x2 h))
           (x4 (+ x3 h)))
      (declare (type real x2 x3 x4))
      (* (/ (* 2 h) 45)
         (+ (*  7 (funcall f x1))
            (* 32 (funcall f x2))
            (* 12 (funcall f x3))
            (* 32 (funcall f x4))
            (*  7 (funcall f x5)))))))

(defun integrate-by-boole (f a b n)
  "Computes and returns the numerical integral of the function F in the
   closed interval [A, B], divided into N equidistant subintervals, and
   obtained by applying Boole's rule to each of the N sections."
  (declare (type (function (real) real) f))
  (declare (type real                   a))
  (declare (type real                   b))
  (declare (type (integer 1 *)          n))
  (let ((k (/ (- b a) n)))
    (declare (type real k))
    (loop
      repeat n
      for x1 of-type real = a then x5
      for x5 of-type real = (+ x1 k)
      sum (calculate-boole-area f x1 x5))))

;;; -------------------------------------------------------

(float
  (integrate-by-boole
    #'(lambda (x) (/ 1 x))
    1 2 8))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; -- Composite Boole's rule.                                      -- ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun integrate-by-composite-booles-rule (f a b n)
  "Computes and returns the numerical integral of the function F in the
   closed interval [A, B], divided into N equidistant subintervals, and
   obtained by following the composite Boole's rule."
  (declare (type (function (real) real) f))
  (declare (type real                   a))
  (declare (type real                   b))
  (declare (type (integer 1 *)          n))
  (let ((h (/ (- b a) n)))
    (declare (type real h))
    (flet ((f[i] (i)
            (declare (type (integer 0 *) i))
            (let ((xi (+ a (* i h))))
              (declare (type real xi))
              (the real (funcall f xi)))))
      (* (/ (* 2 h) 45)
         (+ (*  7 (+ (f[i] 0) (f[i] n)))                                ;; Computes:  7*(y0 + yN).
            (* 32 (loop for i from 1 below n by 2 sum (f[i] i)))        ;; Computes: 32*(y1 + y3 + y5, ...).
            (* 12 (loop for i from 2 below n by 4 sum (f[i] i)))        ;; Computes: 12*(y2 + y6 + y10, ...).
            (* 14 (loop for i from 4 below n by 4 sum (f[i] i))))))))   ;; Computes: 14*(y4 + y8 + y12, ...).

;;; -------------------------------------------------------

(float
  (integrate-by-composite-booles-rule
    #'(lambda (x) (/ 1 x))
    1 2 8))


