/**
 * @file preceiver.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for receiving SCHP5 Packets
 * @version 0.1
 * @date 2025-11-04
 */

 # ifndef PRECEIVER_H
 # define PRECEIVER_H

#include <cstdint>
#include "packet.h"
#include "simple_packet.h"
#include "packet_former.h"
#include "packet_reformer.h"
#include "data_extractor.h"

// State definitions
enum RState {
    RSTART,
    WAIT_DATA,
    CONTROL,
    UNPACKET,
    OUT_BUFF,
    RENDING,
    REND
};

/** 
 * @brief Handles control flow and packeting for SCHP5 data sink
 */
class PacketReceiver {

    private:

    /**
     * @brief The current state of the source control
     */
    RState state;

    int8_t senderAddress;

    uint8_t receiverAddress;

    uint8_t dataLength;

    PacketFormer pFormer;

    PacketReformer pReformer;

    DataExtractor dExtractor;

    float centerFreq;

    uint8_t modulationType;

    uint8_t sequenceNumber;

    /**
     * @brief The buffer of data (from the receiver) to process into packets
     */
    std::vector<bool> * inDataBuffer;

    /**
     * @brief A buffer to store control packets (from the sender)
     */
    std::vector<SimpPacket> * controlPacketBuffer;

    /**
     * @brief A buffer to store the current transaction's data packets (from the sender)
     */
    std::vector<SimpPacket> * dataPacketBuffer;

    /**
     * @brief Store the packets for the current transaction
     */
    std::vector<Packet> * dataPackets;

    /**
     * @brief A buffer to store which packets must be retransmitted
     */
    std::vector<uint8_t> * arqBuffer;

    /**
     * @brief A buffer to store which packets are good
     */
    std::vector<uint8_t> * ackBuffer;

    /**
     * @brief A buffer to store data received
     */
    std::vector<uint8_t> * receiveData;

    /**
     * @brief A buffer to store packets to be sent
     */
    std::vector<Packet> * toSend;

    public:

    /**
     * @brief Construct a new packet transmitter
     */
    PacketReceiver(
        uint8_t sender_address,
        uint8_t receiver_address,
        uint8_t data_length,
        float center_freq,
        uint8_t modulation_type,
        std::vector<bool> * in_data_buffer,
        std::vector<SimpPacket> * control_packet_buffer,
        std::vector<SimpPacket> * data_packet_buffer,
        std::vector<uint8_t> * arq_buffer,
        std::vector<uint8_t> * ack_buffer,
        std::vector<uint8_t> * receive_data,
        std::vector<Packet> * to_send);

    /**
     * @brief Do next operation and if necessary update state
     */
    void tick();

    /**
     * @brief Report current state
     * @returns the current state
     */
    RState getState();

    /**
     * @brief Reports the current center frequency
     * @returns The current center frequency
     */
    float getCenterFrequency();

    /**
     * @brief Reports the current modulation type
     * @returns The current modulation type index, per SCHP5
     */
    uint8_t getModulationType();

};
 
 # endif