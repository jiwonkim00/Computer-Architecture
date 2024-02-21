	.option nopic
	.attribute arch, "rv32i2p1_m2p0"
	.attribute unaligned_access, 0
	.attribute stack_align, 16
	.text
	.align	2
	.globl	digitsum
	.type	digitsum, @function
	
digitsum:
	add a3, a0, x0
	addi a0, x0, 0
	add t0, a0, x0
	addi t1, x0, 10
	
lhs_loop:
	beq a1, x0, rhs_loop
	rem a2, a1, t1
	add t0, t0, a2
	div a1, a1, t1
	jal x0, lhs_loop
	
rhs_loop:
	beq a3, x0, done
	rem a2, a3, t1
	add t0, t0, a2
	div a3, a3, t1
	jal x0, rhs_loop
	
done:
	add a0, t0, x0
	

	#------Your code starts here------
	#LHS: a0, RHS: a1
	
	

	#Load return value to reg a0
	#------Your code ends here------

	#Ret
	jr	ra
	.size	digitsum, .-digitsum
	.ident	"GCC: (g2ee5e430018) 12.2.0"
