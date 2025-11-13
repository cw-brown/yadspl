/**
 * @file ptransmitter.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for transmitting SCHP5 Packets
 * @version 0.1
 * @date 2025-11-04
 */
# ifndef PTRANSMITTER_H
# define PTRANSMITTER_H

#include <cstdint>
#include "packet.h"
#include "simple_packet.h"
#include "packet_former.h"
#include "packet_reformer.h"
#include "data_extractor.h"

// State definitions
enum TState {
    START,
    WAIT_ST,
    WAIT_DATA,
    FORM,
    SEND,
    WAIT_ACK,
    ENDING,
    END
};

/** 
 * @brief Handles control flow and packeting for SCHP5 data source
 */
class PacketTransmitter {

    private:

    /**
     * @brief The current state of the source control
     */
    TState state;

    int8_t senderAddress;

    uint8_t receiverAddress;

    uint8_t dataLength;

    PacketFormer pFormer;

    PacketReformer pReformer;

    float centerFreq;

    uint8_t modulationType;

    uint8_t sequenceNumber;

    Packet nextSend;

    /**
     * @brief The buffer of data (from the receiver) to process into packets
     */
    std::vector<bool> * inDataBuffer;

    /**
     * @brief A buffer to store control packets (from the receiver)
     */
    std::vector<SimpPacket> * controlPacketBuffer;

    /**
     * @brief A buffer to store the current transaction's data packets (from the receiver)
     */
    std::vector<SimpPacket> * dataPacketBuffer;

    /**
     * @brief Store the packets for the current transaction
     */
    std::vector<Packet> * dataPackets;

    /**
     * @brief A buffer to store data to be sent (input)
     */
    std::vector<uint8_t> * toSendData;

    /**
     * @brief A buffer to store packets to be sent (output)
     */
    std::vector<Packet> * toSend;

    public:

    /**
     * @brief Construct a new packet transmitter
     */
    PacketTransmitter(
        uint8_t sender_address,
        uint8_t receiver_address,
        uint8_t data_length,
        float center_freq,
        uint8_t modulation_type,
        std::vector<bool> * in_data_buffer,
        std::vector<SimpPacket> * control_packet_buffer,
        std::vector<SimpPacket> * data_packet_buffer,
        std::vector<uint8_t> * to_send_data,
        std::vector<Packet> * to_send);

    /**
     * @brief Do next operation and if necessary update state
     */
    void tick();

    /**
     * @brief Immediately set the state to ending, to be processed on the next tick
     */
    void end();

    /**
     * @brief Report current state
     * @returns the current state
     */
    TState getState();

};

# endif