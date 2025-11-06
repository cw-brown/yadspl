#include "ptransmitter.h"

PacketTransmitter::PacketTransmitter(
    uint8_t sender_address,
    uint8_t receiver_address,
    uint8_t data_length,
    float center_freq,
    uint8_t modulation_type,
    std::vector<bool> * in_data_buffer,
    std::vector<SimpPacket> * control_packet_buffer,
    std::vector<SimpPacket> * data_packet_buffer,
    std::vector<uint8_t> * to_send_data,
    std::vector<Packet> * to_send) :
    senderAddress(sender_address),
    receiverAddress(receiver_address),
    dataLength(data_length),
    pFormer(sender_address, receiver_address),
    pReformer(in_data_buffer, control_packet_buffer, data_packet_buffer),
    toSendData(to_send_data),
    toSend(to_send),
    state(START),
    sequenceNumber(0) {

    dataPackets = new std::vector<Packet>;
    dataPackets->resize(sizeof(uint8_t));

}

void PacketTransmitter::tick() {

    // Always process whatever has been received, also do whatever is appropriate for state
    pReformer.formPacket(true, senderAddress, receiverAddress);

    switch (state) {

        case START:

            toSend->push_back(pFormer.formTransactionStartPacket(centerFreq, modulationType));
            state = WAIT_ST;
            break;
        
        case WAIT_ST:

            if(!controlPacketBuffer->empty() && (*(controlPacketBuffer->end()-1)).controlCode == 0b1000) {

                state = WAIT_DATA;

            }
            break;

        case WAIT_DATA:

            if(toSendData->size() >= dataLength) {

                state = FORM;

            }
            break;

        case FORM:
            
            dataPackets->at(sequenceNumber) = (pFormer.formDataPacket(toSendData, sequenceNumber, dataLength));
            nextSend = (dataPackets->at(sequenceNumber));
            state = SEND;
            break;

        case SEND:
            
            toSend->push_back(nextSend);
            state = WAIT_ACK;
            break;

        case WAIT_ACK:
            
            // Move on if ack, resend if arq, wait if nothing
            if(!controlPacketBuffer->empty() && (*(controlPacketBuffer->end()-1)).controlCode == 0b0000) {

                // Ack received, move on
                sequenceNumber = (*(controlPacketBuffer->end()-1)).sequenceNumber + 1;
                controlPacketBuffer->pop_back();

                // If sequence if 0, we have rolled over, need to send transaction restart
                if(sequenceNumber == 0) {

                    toSend->push_back(pFormer.formTransactionRestartPacket(centerFreq, modulationType));
                    state = WAIT_ST;

                } else {

                    state = WAIT_DATA;
                    
                }

            } else if(!controlPacketBuffer->empty() && (*(controlPacketBuffer->end()-1)).controlCode == 0b0001) {

                // Arq received, resend indicated packet
                nextSend = dataPackets->at((*(controlPacketBuffer->end()-1)).sequenceNumber);
                controlPacketBuffer->pop_back();
                state = SEND;

            }
            break;

        default:
            
            break;

    }

}