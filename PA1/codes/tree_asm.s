	.option nopic
	.attribute arch, "rv32i2p1_m2p0"
	.attribute unaligned_access, 0
	.attribute stack_align, 16
	.text
	.align	2
	.globl	tree
	.type	tree, @function

tree:

    # Initialize sum to 0
    addi  t0, x0, 0             # t0 = sum

    # Initialize head to the address of the queue (a1)
    addi  t2, a1, 0     # t2 = head

    # Initialize cur to the address of the first element after the queue (a1 + TREE_NODE_SIZE)
    addi  t1, a1, 12  # t1 = cur

    # Copy the root node to the queue
    lw    t3, 0(a0)  # Load root->val
    sw    t3, 0(t2)  # Store val in queue[0]
    lw    t4, 4(a0) # Load root->left
    sw    t4, 4(t2) # Store left in queue[cur]
    lw    t5, 8(a0) # Load root->right
    sw    t5, 8(t2) # Store right in queue[cur]
    

loop:
    # Check if head != cur
    beq   t2, t1, done       # If head == cur, exit the loop
    
    # Load the node at head from the queue
    lw    a1, 0(t2) # Load node->val
    add   t0, t0, a1        # sum += node->val

    # Load node's left child
    lw    a0, 4(t2) # a0 = pointer of left child

    # Check if left child is not NULL
    beq   a0, x0, skip_left      # If node->left is NULL, skip the left child
    
    lw    t3, 0(a0)          # t4 : Get value of the left child (a0: pointer of left child)
    sw    t3, 0(t1)          # Store value of left child in queue[cur]
    lw    t4, 4(a0)          # t4 : Get value of the left child (a0: pointer of left child)
    sw    t4, 4(t1)          # Store value of left child in queue[cur]
    lw    t5, 8(a0)          # t4 : Get value of the left child (a0: pointer of left child)
    sw    t5, 8(t1)          # Store value of left child in queue[cur]
    addi  t1, t1, 12  	     # Increment cur by TREE_NODE_SIZE

skip_left:
    # Load node's right child
    lw    a0, 8(t2) # a0 = pointer of right child

    # Check if right child is not NULL
    beq   a0, x0, skip_right     # If node->right is NULL, skip the right child
    lw    t3, 0(a0)          # t4 : Get value of the left child (a0: pointer of left child)
    sw    t3, 0(t1)          # Store value of left child in queue[cur]
    lw    t4, 4(a0)          # t4 : Get value of the left child (a0: pointer of left child)
    sw    t4, 4(t1)          # Store value of left child in queue[cur]
    lw    t5, 8(a0)          # t4 : Get value of the left child (a0: pointer of left child)
    sw    t5, 8(t1)          # Store value of left child in queue[cur]
    addi  t1, t1, 12  # Increment cur by TREE_NODE_SIZE

skip_right:
    addi  t2, t2, 12  # Increment head by TREE_NODE_SIZE
    jal   x0, loop               # Repeat the loop

done:
    # Move the final result from t0 to a0
    add   a0, t0, x0
    
    
    jr ra
    

	.size	tree, .-tree
	.ident	"GCC: (g2ee5e430018) 12.2.0"
