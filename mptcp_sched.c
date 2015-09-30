/* MPTCP Scheduler module selector. Highly inspired by tcp_cong.c */

#include <linux/module.h>
#include <net/mptcp.h>

static DEFINE_SPINLOCK(mptcp_sched_list_lock);
static LIST_HEAD(mptcp_sched_list);

struct defsched_priv {
	u32	last_rbuf_opti;
};

static struct defsched_priv *defsched_get_priv(const struct tcp_sock *tp)
{
	return (struct defsched_priv *)&tp->mptcp->mptcp_sched[0];
}

/* If the sub-socket sk available to send the skb? */
static bool mptcp_is_available(struct sock *sk, struct sk_buff *skb,
			       bool zero_wnd_test)
{
	struct tcp_sock *tp = tcp_sk(sk);
	unsigned int mss_now, space, in_flight;

	/* Set of states for which we are allowed to send data */
	if (!mptcp_sk_can_send(sk))
		return false;

	/* We do not send data on this subflow unless it is
	 * fully established, i.e. the 4th ack has been received.
	 */
	if (tp->mptcp->pre_established)
		return false;

	if (tp->pf)
		return false;

	if (inet_csk(sk)->icsk_ca_state == TCP_CA_Loss) {
		/* If SACK is disabled, and we got a loss, TCP does not exit
		 * the loss-state until something above high_seq has been acked.
		 * (see tcp_try_undo_recovery)
		 *
		 * high_seq is the snd_nxt at the moment of the RTO. As soon
		 * as we have an RTO, we won't push data on the subflow.
		 * Thus, snd_una can never go beyond high_seq.
		 */
		if (!tcp_is_reno(tp))
			return false;
		else if (tp->snd_una != tp->high_seq)
			return false;
	}

	if (!tp->mptcp->fully_established) {
		/* Make sure that we send in-order data */
		if (skb && tp->mptcp->second_packet &&
		    tp->mptcp->last_end_data_seq != TCP_SKB_CB(skb)->seq)
			return false;
	}

	/* If TSQ is already throttling us, do not send on this subflow. When
	 * TSQ gets cleared the subflow becomes eligible again.
	 */
	if (test_bit(TSQ_THROTTLED, &tp->tsq_flags))
		return false;

	in_flight = tcp_packets_in_flight(tp);
	/* Not even a single spot in the cwnd */
	if (in_flight >= tp->snd_cwnd)
		return false;

	/* Now, check if what is queued in the subflow's send-queue
	 * already fills the cwnd.
	 */
	space = (tp->snd_cwnd - in_flight) * tp->mss_cache;

	if (tp->write_seq - tp->snd_nxt > space)
		return false;

	if (zero_wnd_test && !before(tp->write_seq, tcp_wnd_end(tp)))
		return false;

	mss_now = tcp_current_mss(sk);

	/* Don't send on this subflow if we bypass the allowed send-window at
	 * the per-subflow level. Similar to tcp_snd_wnd_test, but manually
	 * calculated end_seq (because here at this point end_seq is still at
	 * the meta-level).
	 */
	if (skb && !zero_wnd_test &&
	    after(tp->write_seq + min(skb->len, mss_now), tcp_wnd_end(tp)))
		return false;

	return true;
}

/* Are we not allowed to reinject this skb on tp? */
static int mptcp_dont_reinject_skb(struct tcp_sock *tp, struct sk_buff *skb)
{
	/* If the skb has already been enqueued in this sk, try to find
	 * another one.
	 */
	return skb &&
		/* Has the skb already been enqueued into this subsocket? */
		mptcp_pi_to_flag(tp->mptcp->path_index) & TCP_SKB_CB(skb)->path_mask;
}


/*
added by bh oh
this fuction is used for CP scheduling and 
it collects the subflow information in MPTCP layer.
*/
int CP_update_subflow_info(struct CP_subflow_info *cp_info_p, struct mptcp_cb *mpcb)				   
{

	struct sock *sk;
	int subflow_N=0;
	int i=0,j=0;

	mptcp_for_each_sk(mpcb, sk) {
		struct tcp_sock *tp1 = tcp_sk(sk);
		
	// this is initialization of CP scheduling
	// the reason that the initialization is located in this function is 
	// due to the SIGSEGV error when it operated in ns-3-dce-3.14  
		if (tp1->initial_CP != 100 )
		{  
			tp1->blocking_count = 0;
			tp1->prevention_time_by_FP = 0;
			tp1->blocking_RTT = 0;
			tp1->blocking_seq = 0;
			tp1->prevent_by_FP = 0;
			tp1->total_T_guard = 0;
			tp1->prob_count = 0;
			tp1->initial_CP = 100;
		}

// Note that below annotated codes are only for the NS-3-dce examples (like dce-cradle-mptcp, 
// and dce-iperf-mptcp) in ns-3-dce-3.11.
// becuase in these example (like dce-cradle-mptcp, and dce-iperf-mptcp), 
// inappropriate routing can be generated.
// although inappropriate routing does not causes the serious problems, 
// this code may correct the inappropriate transmission. by bh oh.
///*
		if ( tp1->mptcp->path_index == 2 || tp1->mptcp->path_index == 3 )
		{
			cp_info_p->RTT_subflow[subflow_N] = 9999;
		} 
		else
		{
			if ( tp1->srtt == 0 )
			{
				cp_info_p->RTT_subflow[subflow_N] = 9999;
			}
			else
			{
				cp_info_p->RTT_subflow[subflow_N] = tp1->srtt;
			}
		}
//*/ 
// this codes should be disabled when above codes are enalbed.
/*
		if ( tp1->srtt == 0 )
		{
			cp_info_p->RTT_subflow[subflow_N] = 9999;
		}
		else
		{
			cp_info_p->RTT_subflow[subflow_N] = tp1->srtt;
		}
*/

		cp_info_p->RTT_subflow_index[subflow_N] = tp1->mptcp->path_index;
		cp_info_p->CWND_subflow[subflow_N] = tp1->snd_cwnd;
		cp_info_p->outstanding_subflow[subflow_N] = tp1->packets_out + tp1->retrans_out;
		cp_info_p->subflow_mss[subflow_N] = tp1->advmss;
		cp_info_p->subflow_ssthresh[subflow_N] = tp1->snd_ssthresh;

		subflow_N=subflow_N+1;

	}
	// sorting the subflow info
	for (i=0; i<subflow_N-1; i++)
	{
		for (j=i; j<subflow_N-1; j++)
		{
			if ( cp_info_p->RTT_subflow[i] > cp_info_p->RTT_subflow[j+1] )
			{
				u64 RTT_box;
				u64 index_box;
				u64 cwnd_box;				
				u64 outstanding_box;
				u64 mss_box;
				u64 mss_ssthresh;				
				
				RTT_box = cp_info_p->RTT_subflow[i];
				cp_info_p->RTT_subflow[i] = cp_info_p->RTT_subflow[j+1];
				cp_info_p->RTT_subflow[j+1] = RTT_box;

				index_box = cp_info_p->RTT_subflow_index[i];
				cp_info_p->RTT_subflow_index[i] = cp_info_p->RTT_subflow_index[j+1];
				cp_info_p->RTT_subflow_index[j+1] = index_box;

				cwnd_box = cp_info_p->CWND_subflow[i];
				cp_info_p->CWND_subflow[i] = cp_info_p->CWND_subflow[j+1];
				cp_info_p->CWND_subflow[j+1] = cwnd_box;

				outstanding_box = cp_info_p->outstanding_subflow[i];
				cp_info_p->outstanding_subflow[i] = cp_info_p->outstanding_subflow[j+1];
				cp_info_p->outstanding_subflow[j+1] = outstanding_box;

				mss_box = cp_info_p->subflow_mss[i];
				cp_info_p->subflow_mss[i]= cp_info_p->subflow_mss[j+1];
				cp_info_p->subflow_mss[j+1] = mss_box;
				
				mss_ssthresh = cp_info_p->subflow_ssthresh[i];
				cp_info_p->subflow_ssthresh[i] = cp_info_p->subflow_ssthresh[j+1];
				cp_info_p->subflow_ssthresh[j+1] = mss_ssthresh;
			}
		}
	}



	return subflow_N;

}



/* This is the scheduler. This function decides on which flow to send
 * a given MSS. If all subflows are found to be busy, NULL is returned
 * The flow is selected based on the shortest RTT.
 * If all paths have full cong windows, we simply return NULL.
 *
 * Additionally, this function is aware of the backup-subflows.
 */
static struct sock *get_available_subflow(struct sock *meta_sk,
					  struct sk_buff *skb,
					  bool zero_wnd_test)
{
	struct mptcp_cb *mpcb = tcp_sk(meta_sk)->mpcb;
	struct sock *sk, *bestsk = NULL, *lowpriosk = NULL, *backupsk = NULL;
	u32 min_time_to_peer = 0xffffffff, lowprio_min_time_to_peer = 0xffffffff;
	int cnt_backups = 0;
	int subflow_N=0;
	int i=0, j=0, k=0, l=0, m=0;
	int path_change=0;

	struct tcp_sock *meta_tp = tcp_sk(meta_sk);

	/* if there is only one subflow, bypass the scheduling function */
	if (mpcb->cnt_subflows == 1) {
		bestsk = (struct sock *)mpcb->connection_list;
		if (!mptcp_is_available(bestsk, skb, zero_wnd_test))
			bestsk = NULL;
		return bestsk;
	}

	/* Answer data_fin on same subflow!!! */
	if (meta_sk->sk_shutdown & RCV_SHUTDOWN &&
	    skb && mptcp_is_data_fin(skb)) {
		mptcp_for_each_sk(mpcb, sk) {
			if (tcp_sk(sk)->mptcp->path_index == mpcb->dfin_path_index &&
			    mptcp_is_available(sk, skb, zero_wnd_test))
				return sk;
		}
	}

	// added by bh oh
	struct CP_subflow_info cp_info;
	struct CP_subflow_info *cp_info_p = &cp_info;
	subflow_N = CP_update_subflow_info(cp_info_p, mpcb); 
	// added end

	/* First, find the best subflow */
	mptcp_for_each_sk(mpcb, sk) {
		struct tcp_sock *tp = tcp_sk(sk);
		int this_mss;
		int this_mss1;
		

// Note that below annotated codes are only for the NS-3-dce examples (like dce-cradle-mptcp, 
// and dce-iperf-mptcp) in ns-3-dce-3.14.
// becuase in these example (like dce-cradle-mptcp, and dce-iperf-mptcp), 
// inappropriate routing can be generated.
// although inappropriate routing does not causes the serious problems, 
// this code may correct the inappropriate transmission. by bh oh.
/*   //----------------------------------------------------------------------
		if ( tp->mptcp->path_index == 2 || tp->mptcp->path_index == 3   )
		{
			continue;
		}
*/  //------------------------------------------------------------------------ 

// /* annotation point for disabling CP-scheduling


		if ( tp->mptcp->path_index == cp_info_p->RTT_subflow_index[0] )
		{ // this path is fastest path, just pass

			// for conservation scheudling
/*
			if (cp_info_p->outstanding_subflow[1]*cp_info_p->subflow_mss[1]*2 >= (meta_tp->snd_wnd - cp_info_p->subflow_mss[1]))
			{
				printk("just uses previous fastest path. \n");
				continue;
			}
*/
		}

// for conservation scheudling
/*
		else if ( (cp_info_p->CWND_subflow[0] < cp_info_p->subflow_ssthresh[0]) || (cp_info_p->outstanding_subflow[0]*cp_info_p->subflow_mss[0]*2 >= (meta_tp->snd_wnd - cp_info_p->subflow_mss[0])) )
		{	// for use only fastest path unless sufficient buffer is available
			int sk_index=0;
			for (k=0; k<subflow_N; k++)
			{
				if ( tp->mptcp->path_index == cp_info_p->RTT_subflow_index[k] )
				{
					sk_index=k;
				}
			}
			if ( cp_info_p->RTT_subflow[sk_index] != 9999 )			
			{
				tp->prevent_by_FP=0; // stop the probing
				continue;		
			}
		}
*/
		else
		{
			int sk_index=0;
			for (k=0; k<subflow_N; k++)
			{
				if ( tp->mptcp->path_index == cp_info_p->RTT_subflow_index[k] )
				{
					sk_index=k; // find the RTT rank in this path.
				}
			}

			u64 estimate_out=0;
			u64 estimate_out_byte=0;
			u64 estimate_out_slow=0;
			u64 estimate_out_byte_slow=0;
			u64 CWND_slow=0;
			u64 sk_RTT;

			u32 rto = inet_csk(sk)->icsk_rto;

			if ( cp_info_p->RTT_subflow[sk_index] != 9999 )
			{
				for(l=0; l<sk_index; l++)
				{
					estimate_out=0;
					sk_RTT = cp_info_p->RTT_subflow[sk_index];
					int count=0;
//					printk("sk_RTT is  %ld \n", sk_RTT);
//					printk("RTT_subflow[l] is  %ld \n", RTT_subflow[l]);

					while( sk_RTT > cp_info_p->RTT_subflow[l] )
					{
						count=count+1;
						sk_RTT= sk_RTT - cp_info_p->RTT_subflow[l];
					}
					for (m=1; m<count; m++) 
					{
						estimate_out=estimate_out + cp_info_p->outstanding_subflow[l];				
					}
					estimate_out_byte= estimate_out_byte +(estimate_out*cp_info_p->subflow_mss[l]);
					estimate_out_slow = cp_info_p->outstanding_subflow[sk_index];
					estimate_out_byte_slow = estimate_out_slow*cp_info_p->subflow_mss[sk_index];
					CWND_slow = cp_info_p->CWND_subflow[sk_index]*cp_info_p->subflow_mss[sk_index];
				} // calculate the required buffer size to handle the out-of-order pkts

				if (skb)
				{
					if ( mptcp_is_available(sk, skb, zero_wnd_test) )
					{
						u32 end_seq1 = TCP_SKB_CB(skb)->end_seq;
						u32 meta_tp_wnd_end = meta_tp->snd_una + meta_tp->snd_wnd;
						u32 available_buffer_byte;
						if ( meta_tp_wnd_end > end_seq1 )
						{
							available_buffer_byte=meta_tp_wnd_end - end_seq1;
						}
						else
						{
							available_buffer_byte=0;
						}

						if ( estimate_out_byte > available_buffer_byte )
						{	
								u32 rtt_time_for_blocking;
								int beta=4; // for initial path probing period
								if (tp->blocking_count == 0 )
								{
									rtt_time_for_blocking=(tp->srtt >> 3)*beta; 
								}
								else
								{
									rtt_time_for_blocking=tp->blocking_RTT*beta;
	//								rtt_time_for_blocking=rtt_time_for_blocking*tp->blocking_count; // for dynamic probing
								}	

								
								if ( rtt_time_for_blocking > 2500 )
								{ // this is maximum path probing period
									rtt_time_for_blocking = 2500;
									printk("probing duration limitation! \n");
								}


								if ( (tcp_time_stamp > tp->prevention_time_by_FP+ rtt_time_for_blocking) && (tp->blocking_seq < meta_tp->snd_una) )
								{
									tp->prevent_by_FP=10; // the amount of path probings. 
									tp->prevention_time_by_FP = tcp_time_stamp;
								
									tp->blocking_seq=meta_tp->snd_una;
									tp->blocking_RTT=(tp->srtt >> 3);
									tp->blocking_count=tp->blocking_count+1;
									tp->total_T_guard+=rtt_time_for_blocking;
	//								printk("tp->blocking_count is %ld \n", tp->blocking_count);
	//								printk("tp->total_T_guard is %ld \n", tp->total_T_guard);
								
									continue;
								}
								else if ( tp->blocking_count > 1)
								{
									if ( tp->prevent_by_FP <= 0 )
									{
										continue;					
									}
									else
									{
										tp->prevent_by_FP=tp->prevent_by_FP-1;
										tp->prob_count+=1;
										// operate the path probing
									}
								}
								else
								{
									continue;
								}
						}

						else // when buffer size > estimate_out_byte
						{
							u64 available_buffer=0;
							u64 total_outstanding_byte=0;

							for (k=0; k<subflow_N; k++)
							{
								total_outstanding_byte+= cp_info_p->outstanding_subflow[k]*cp_info_p->subflow_mss[k];
							}

							available_buffer = available_buffer_byte - estimate_out_byte;

							if ( available_buffer < total_outstanding_byte )
							{
//								continue;  // for conservative scheduling
							}
							/*
							else if ( cp_info_p->RTT_subflow[sk_index] > 2*cp_info_p->RTT_subflow[0] ) 
							{ // for delay constraints			
								meta_tp->blocking_for_delay+=1;
								continue;  // for conservative scheduling with delay constraint
							}
							*/
						}
					}
				}
			}				
		}




// */ annotation point for disabling CP-scheduling



		if (tp->mptcp->rcv_low_prio || tp->mptcp->low_prio)
			cnt_backups++;

		if ((tp->mptcp->rcv_low_prio || tp->mptcp->low_prio) &&
		    tp->srtt < lowprio_min_time_to_peer) {
			if (!mptcp_is_available(sk, skb, zero_wnd_test))
				continue;

			if (mptcp_dont_reinject_skb(tp, skb)) {
				backupsk = sk;
				continue;
			}

			lowprio_min_time_to_peer = tp->srtt;
			lowpriosk = sk;
		} else if (!(tp->mptcp->rcv_low_prio || tp->mptcp->low_prio) &&
			   tp->srtt < min_time_to_peer) {
			if (!mptcp_is_available(sk, skb, zero_wnd_test))
				continue;

			if (mptcp_dont_reinject_skb(tp, skb)) {
				backupsk = sk;
				continue;
			}

			min_time_to_peer = tp->srtt;
			bestsk = sk;
		}
	}

	if (mpcb->cnt_established == cnt_backups && lowpriosk) {
		sk = lowpriosk;
	} else if (bestsk) {
		sk = bestsk;
	} else if (backupsk) {
		/* It has been sent on all subflows once - let's give it a
		 * chance again by restarting its pathmask.
		 */
		if (skb)
			TCP_SKB_CB(skb)->path_mask = 0;
		sk = backupsk;
	}

	return sk;
}

static struct sk_buff *mptcp_rcv_buf_optimization(struct sock *sk, int penal)
{
	struct sock *meta_sk;
	struct tcp_sock *tp = tcp_sk(sk), *tp_it;
	struct sk_buff *skb_head;
	struct defsched_priv *dsp = defsched_get_priv(tp);

	if (tp->mpcb->cnt_subflows == 1)
		return NULL;

	meta_sk = mptcp_meta_sk(sk);
	skb_head = tcp_write_queue_head(meta_sk);

	if (!skb_head || skb_head == tcp_send_head(meta_sk))
		return NULL;

	/* If penalization is optional (coming from mptcp_next_segment() and
	 * We are not send-buffer-limited we do not penalize. The retransmission
	 * is just an optimization to fix the idle-time due to the delay before
	 * we wake up the application.
	 */
	if (!penal && sk_stream_memory_free(meta_sk))
		goto retrans;

	/* Only penalize again after an RTT has elapsed */
	if (tcp_time_stamp - dsp->last_rbuf_opti < tp->srtt >> 3)
		goto retrans;

	/* Half the cwnd of the slow flow */
	mptcp_for_each_tp(tp->mpcb, tp_it) {
		if (tp_it != tp &&
		    TCP_SKB_CB(skb_head)->path_mask & mptcp_pi_to_flag(tp_it->mptcp->path_index)) {
			if (tp->srtt < tp_it->srtt && inet_csk((struct sock *)tp_it)->icsk_ca_state == TCP_CA_Open) {
				u32 prior_cwnd = tp_it->snd_cwnd;

				tp_it->snd_cwnd = max(tp_it->snd_cwnd >> 1U, 1U);

				/* If in slow start, do not reduce the ssthresh */
				if (prior_cwnd >= tp_it->snd_ssthresh)
					tp_it->snd_ssthresh = max(tp_it->snd_ssthresh >> 1U, 2U);

				dsp->last_rbuf_opti = tcp_time_stamp;
			}
			break;
		}
	}

retrans:

	/* Segment not yet injected into this path? Take it!!! */
	if (!(TCP_SKB_CB(skb_head)->path_mask & mptcp_pi_to_flag(tp->mptcp->path_index))) {
		bool do_retrans = false;
		mptcp_for_each_tp(tp->mpcb, tp_it) {
			if (tp_it != tp &&
			    TCP_SKB_CB(skb_head)->path_mask & mptcp_pi_to_flag(tp_it->mptcp->path_index)) {
				if (tp_it->snd_cwnd <= 4) {
					do_retrans = true;
					break;
				}

				if (4 * tp->srtt >= tp_it->srtt) {
					do_retrans = false;
					break;
				} else {
					do_retrans = true;
				}
			}
		}

		if (do_retrans && mptcp_is_available(sk, skb_head, false))
			return skb_head;
	}
	return NULL;
}

/* Returns the next segment to be sent from the mptcp meta-queue.
 * (chooses the reinject queue if any segment is waiting in it, otherwise,
 * chooses the normal write queue).
 * Sets *@reinject to 1 if the returned segment comes from the
 * reinject queue. Sets it to 0 if it is the regular send-head of the meta-sk,
 * and sets it to -1 if it is a meta-level retransmission to optimize the
 * receive-buffer.
 */
static struct sk_buff *__mptcp_next_segment(struct sock *meta_sk, int *reinject)
{
	struct mptcp_cb *mpcb = tcp_sk(meta_sk)->mpcb;
	struct sk_buff *skb = NULL;

	*reinject = 0;

	/* If we are in fallback-mode, just take from the meta-send-queue */
	if (mpcb->infinite_mapping_snd || mpcb->send_infinite_mapping)
		return tcp_send_head(meta_sk);

	skb = skb_peek(&mpcb->reinject_queue);

	if (skb) {
		*reinject = 1;
	} else {
		skb = tcp_send_head(meta_sk);

		if (!skb && meta_sk->sk_socket &&
		    test_bit(SOCK_NOSPACE, &meta_sk->sk_socket->flags) &&
		    sk_stream_wspace(meta_sk) < sk_stream_min_wspace(meta_sk)) {
			struct sock *subsk = get_available_subflow(meta_sk, NULL,
								   false);
			if (!subsk)
				return NULL;

//			skb = mptcp_rcv_buf_optimization(subsk, 0);
			skb = NULL; // modified by bh oh
			if (skb)
				*reinject = -1;
		}
	}
	return skb;
}

static struct sk_buff *mptcp_next_segment(struct sock *meta_sk,
					  int *reinject,
					  struct sock **subsk,
					  unsigned int *limit)
{
	struct sk_buff *skb = __mptcp_next_segment(meta_sk, reinject);
	unsigned int mss_now;
	struct tcp_sock *subtp;
	u16 gso_max_segs;
	u32 max_len, max_segs, window, needed;

	/* As we set it, we have to reset it as well. */
	*limit = 0;

	if (!skb)
		return NULL;

	*subsk = get_available_subflow(meta_sk, skb, false);
	if (!*subsk)
		return NULL;

	subtp = tcp_sk(*subsk);
	mss_now = tcp_current_mss(*subsk);

	if (!*reinject && unlikely(!tcp_snd_wnd_test(tcp_sk(meta_sk), skb, mss_now))) {
//		skb = mptcp_rcv_buf_optimization(*subsk, 1);
		skb = NULL; // modified by bh oh
		if (skb)
			*reinject = -1;
		else
			return NULL;
	}

	/* No splitting required, as we will only send one single segment */
	if (skb->len <= mss_now)
		return skb;

	/* The following is similar to tcp_mss_split_point, but
	 * we do not care about nagle, because we will anyways
	 * use TCP_NAGLE_PUSH, which overrides this.
	 *
	 * So, we first limit according to the cwnd/gso-size and then according
	 * to the subflow's window.
	 */

	gso_max_segs = (*subsk)->sk_gso_max_segs;
	if (!gso_max_segs) /* No gso supported on the subflow's NIC */
		gso_max_segs = 1;
	max_segs = min_t(unsigned int, tcp_cwnd_test(subtp, skb), gso_max_segs);
	if (!max_segs)
		return NULL;

	max_len = mss_now * max_segs;
	window = tcp_wnd_end(subtp) - subtp->write_seq;

	needed = min(skb->len, window);
	if (max_len <= skb->len)
		/* Take max_win, which is actually the cwnd/gso-size */
		*limit = max_len;
	else
		/* Or, take the window */
		*limit = needed;

	return skb;
}

static void defsched_init(struct sock *sk)
{
	struct defsched_priv *dsp = defsched_get_priv(tcp_sk(sk));

	dsp->last_rbuf_opti = tcp_time_stamp;
}

struct mptcp_sched_ops mptcp_sched_default = {
	.get_subflow = get_available_subflow,
	.next_segment = mptcp_next_segment,
	.init = defsched_init,
	.name = "default",
	.owner = THIS_MODULE,
};

static struct mptcp_sched_ops *mptcp_sched_find(const char *name)
{
	struct mptcp_sched_ops *e;

	list_for_each_entry_rcu(e, &mptcp_sched_list, list) {
		if (strcmp(e->name, name) == 0)
			return e;
	}

	return NULL;
}

int mptcp_register_scheduler(struct mptcp_sched_ops *sched)
{
	int ret = 0;

	if (!sched->get_subflow || !sched->next_segment)
		return -EINVAL;

	spin_lock(&mptcp_sched_list_lock);
	if (mptcp_sched_find(sched->name)) {
		pr_notice("%s already registered\n", sched->name);
		ret = -EEXIST;
	} else {
		list_add_tail_rcu(&sched->list, &mptcp_sched_list);
		pr_info("%s registered\n", sched->name);
	}
	spin_unlock(&mptcp_sched_list_lock);

	return ret;
}
EXPORT_SYMBOL_GPL(mptcp_register_scheduler);

void mptcp_unregister_scheduler(struct mptcp_sched_ops *sched)
{
	spin_lock(&mptcp_sched_list_lock);
	list_del_rcu(&sched->list);
	spin_unlock(&mptcp_sched_list_lock);
}
EXPORT_SYMBOL_GPL(mptcp_unregister_scheduler);

void mptcp_get_default_scheduler(char *name)
{
	struct mptcp_sched_ops *sched;

	BUG_ON(list_empty(&mptcp_sched_list));

	rcu_read_lock();
	sched = list_entry(mptcp_sched_list.next, struct mptcp_sched_ops, list);
	strncpy(name, sched->name, MPTCP_SCHED_NAME_MAX);
	rcu_read_unlock();
}

int mptcp_set_default_scheduler(const char *name)
{
	struct mptcp_sched_ops *sched;
	int ret = -ENOENT;

	spin_lock(&mptcp_sched_list_lock);
	sched = mptcp_sched_find(name);
#ifdef CONFIG_MODULES
	if (!sched && capable(CAP_NET_ADMIN)) {
		spin_unlock(&mptcp_sched_list_lock);

		request_module("mptcp_%s", name);
		spin_lock(&mptcp_sched_list_lock);
		sched = mptcp_sched_find(name);
	}
#endif

	if (sched) {
		list_move(&sched->list, &mptcp_sched_list);
		ret = 0;
	} else {
		pr_info("%s is not available\n", name);
	}
	spin_unlock(&mptcp_sched_list_lock);

	return ret;
}

void mptcp_init_scheduler(struct mptcp_cb *mpcb)
{
	struct mptcp_sched_ops *sched;

	rcu_read_lock();
	list_for_each_entry_rcu(sched, &mptcp_sched_list, list) {
		if (try_module_get(sched->owner)) {
			mpcb->sched_ops = sched;
			break;
		}
	}
	rcu_read_unlock();
}

/* Manage refcounts on socket close. */
void mptcp_cleanup_scheduler(struct mptcp_cb *mpcb)
{
	module_put(mpcb->sched_ops->owner);
}

/* Set default value from kernel configuration at bootup */
static int __init mptcp_scheduler_default(void)
{
	BUILD_BUG_ON(sizeof(struct defsched_priv) > MPTCP_SCHED_SIZE);

	return mptcp_set_default_scheduler(CONFIG_DEFAULT_MPTCP_SCHED);
}
late_initcall(mptcp_scheduler_default);
