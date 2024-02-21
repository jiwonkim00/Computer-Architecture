	.option nopic
	.attribute arch, "rv32i2p1_m2p0"
	.attribute unaligned_access, 0
	.attribute stack_align, 16
	.text
	.align	2
	.globl	fibonacci
	.type	fibonacci, @function
	
fibonacci:
	addi sp, sp, -12
	sw ra, 8(sp)
	sw a0, 0(sp)
	
	addi t0, x0, 2
	lw a1, 0(sp)
	slt t1, t0, a1
	bne t1, x0, not_base_case
	addi a0, x0, 1
	lw ra, 8(sp)
	addi sp, sp, 12
	jalr x0, 0(ra)
	
not_base_case:
	lw a0, 0(sp)
	addi a0, a0, -1
	jal ra, fibonacci
	sw a0, 4(sp)
	
	lw a0, 0(sp)
	
	addi a0, a0, -2
	jal ra, fibonacci
	lw a1, 4(sp)
	add a0, a1, a0
	
	lw ra, 8(sp)
	addi sp, sp, 12
	jalr x0, 0(ra)

	#------Your code starts here------
	# fibonacci number : a0

	# Load return value to reg a0
	#------Your code ends here------

	jr	ra
	.size	fibonacci, .-fibonacci
	.ident	"GCC: (g2ee5e430018) 12.2.0"
